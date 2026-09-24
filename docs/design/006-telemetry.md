# 006: Benchmark Execution Telemetry

[![Status: Implemented](https://img.shields.io/badge/Status-Implemented-green.svg)](https://github.com/omnibenchmark/docs/design)
[![Version: 0.2](https://img.shields.io/badge/Version-0.2-blue.svg)](https://github.com/omnibenchmark/docs/design)

**Authors**: ben
**Date**: 2026-02-19
**Status**: Implemented — [PR #333](https://github.com/omnibenchmark/omnibenchmark/pull/333), merged to `main` as `047172c` (2026-06-23)
**Version**: 0.2
**Supersedes**: N/A
**Reviewed-by**: TBD
**Related Issues**: [#293](https://github.com/omnibenchmark/omnibenchmark/pull/293) (original design PR)

## Changes

| Version | Date | Description | Author |
|---------|------|-------------|--------|
| 0.1 | 2026-02-19 | Initial draft | ben |
| 0.2 | 2026-09-24 | Reviewed against the merged implementation (`047172c`); status Implemented; CLI, relay and status table corrected | ben |

## 1. Problem Statement

Benchmark runs are long-lived, distributed computations. Understanding what happened — which rules ran, how long they took, what failed and why — currently requires reading Snakemake logs, per-rule log files, and stdout. There is no structured, machine-readable record of a benchmark execution that downstream tools (dashboards, CI systems, result databases) can consume.

Additionally, the current `ob run` execution model produces no progress signal compatible with distributed observability infrastructure. Reproducibility and auditability of benchmark runs require richer provenance than what Snakemake natively provides.

## 2. Design Goals

- **Structured execution record**: A benchmark run can produce a machine-readable trace — a hierarchical record of what ran, when, how long, and whether it succeeded.
- **OTLP compatibility**: Output conforms to the OpenTelemetry Protocol (OTLP) so it can be forwarded to any standard tracing backend (Jaeger, Tempo, Aspire Dashboard, etc.) without custom tooling.
- **No infrastructure required**: Output is plain NDJSON on stdout or in a local file.
- **Non-intrusive**: Telemetry is opt-in at the CLI level. Benchmarks that don't use `--telemetry` are unaffected.
- **Explicit-Snakefile dependency is explicit**: Telemetry is only available on the explicit Snakefile execution path (§3.5).

### Non-Goals

- Replacing Snakemake's native logging or `--report`.
- A built-in dashboard. Consumption is delegated to [obmon](https://github.com/omnibenchmark/obmon).
- Automatic upload to a centralised telemetry service.
- Telemetry for non-Snakemake backends.
- Metrics (counters, histograms) — only distributed traces (spans) and logs.

## 3. Specification

### 3.1 Output Format

Telemetry is written as **OTLP JSON Lines** (one JSON object per line, newline-delimited). Each line is either a `ResourceSpans` or `ResourceLogs` envelope as defined by the OTLP/JSON spec.

This format is:
- Self-contained (no schema file needed)
- Appendable (each line is independent)
- Directly ingestible by OTLP-aware backends
- Human-readable with standard JSON tools

Spans are emitted when they **complete**: a rule span when its rule finishes, a module span when all its rules finish, a stage span when all its modules finish, and the benchmark span last. A killed run therefore leaves completed rule spans without a root span.

### 3.2 Span Hierarchy

A benchmark run is represented as a four-level span tree:

```
benchmark: {name}                    # root: one per ob run invocation
├── setup: module resolution         # setup phase spans
├── setup: environment preparation
├── stage: {stage_id}                # one per stage
│   └── module: {module_id}          # one per (stage, module_id)
│       └── rule: {rule_name}        # one per resolved node
│           ├── [event] stdout
│           ├── [event] stderr
│           └── [event] exception    # only on failure
└── ...
```

**Root span**: `benchmark.name`, `benchmark.version`, `benchmark.author`, `benchmark.total_rules`, `benchmark.software_backend`, `benchmark.cores`.

**Stage span**: `stage.id`, `stage.name`, `stage.module_count`.

**Module span**: `module.id`, `module.stage_id`, `module.rule_count`.

**Rule span**:
- `rule.name`, `rule.stage_id`, `rule.module_id`, `rule.node_id`, `rule.param_id`
- `rule.inputs` (list of paths), `rule.outputs` (list of paths)
- `rule.parameters` (JSON-encoded parameter dict)
- `rule.skipped` — set when Snakemake did not re-run the rule (output already up to date)

Rule stdout/stderr is attached twice: as span events and as OTLP log records correlated to the span. Output is captured by piping each rule through `tee` into its `log:` file (`backend/snakemake.py`).

All spans share one `trace_id`, generated per invocation.

**Timing.** Rule start and end are observed by the `ob` process parsing Snakemake's stdout, not measured inside the job. Durations are orchestrator wall clock.

### 3.3 CLI Interface

```
ob run benchmark.yaml [--telemetry] [--telemetry-output PATH]
```

- `--telemetry`: emit telemetry to **stdout**. Disables the Rich progress UI, which would compete for stdout, and suppresses logging.
- `--telemetry-output PATH`: write to `PATH` instead; keeps Rich active. Implies `--telemetry`. The file is truncated on open.

There is no default file location.

### 3.4 Consuming the Stream

The integration contract is the NDJSON itself. [obmon](https://github.com/omnibenchmark/obmon) tails the file (or reads stdout) and provides the dashboard and OTLP forwarding. No relay script ships in this repository, and the core package has no OpenTelemetry SDK dependency.

### 3.5 Dependency on Explicit Snakefile

Telemetry requires the explicit Snakefile execution path. This is the only supported path for `ob run` today, but the dependency is worth documenting:

- **Why:** Telemetry initialises the full span hierarchy (benchmark → stage → module → rule) from the resolved node list *before* Snakemake executes. This is only possible when the complete DAG is materialised upfront (i.e. wildcards resolved, all node IDs known).
- **Consequence:** If a future backend uses native Snakemake wildcard expansion (not explicit rules), telemetry as specified here cannot work without a different integration approach.
- **Constraint captured:** `ob run` always uses `_generate_explicit_snakefile()`. Any backend that bypasses explicit rule generation loses telemetry support.

## 4. Implementation Status

Merged in `047172c` ([PR #333](https://github.com/omnibenchmark/omnibenchmark/pull/333)).

| Component | Status | Notes |
|-----------|--------|-------|
| `TelemetryEmitter` (`telemetry/emitter.py`) | Done | No external deps |
| OTLP JSON output | Done | Spans + correlated log records |
| `--telemetry` / `--telemetry-output` | Done | See §3.3 |
| Span hierarchy (benchmark/stage/module/rule) | Done | Shared `trace_id` |
| Per-rule stdout/stderr capture | Done | Via `tee` into the rule `log:` file |
| Tests | Done | `tests/telemetry/` (emitter, events, spans) |
| External forwarding | Moved out | obmon replaces the proposed `scripts/telemetry-relay.py` |
| `SpanBuilder` (`telemetry/spans.py`) | Unused | Exported and tested, but the emitter does not use it |
| Default output location | Not done | Stdout unless `--telemetry-output` is given |
| Run ID linked to manifest | Not done | `write_run_manifest()` is called without the trace id (`cli/run.py`), so `manifest.json`'s `run_id` differs from `trace_id` despite the docstring |
| Software environment in rule spans | Not done | |
| YAML configuration | Not done | |

## 5. Known Issues and Gaps

**Resume overwrites the trace.** `--telemetry-output` truncates its file. Re-running into the same path after an interruption discards the earlier trace.

**Zero-duration rules.** If a rule's start line is missed in Snakemake's stdout, `rule_completed` sets start = end, producing a 0 s span that looks real. Consumers should drop non-positive durations.

**Cluster timing.** Under a non-local Snakemake executor, the observed "start" is close to job submission, so queue wait is counted as runtime, and the host that ran the job is not recorded.

**Dry runs do not enumerate nodes.** `init_benchmark()` only fills internal state; `ob run --dry --telemetry` emits the module-resolution span and nothing per node.

**Environment metadata missing.** Rule spans don't include which software environment was used, whether it was a container or conda env, or what the resolved image/env path was.

## 6. Planned Improvements

| Item | Description | Priority |
|------|-------------|----------|
| Run ID | Pass the trace id to `write_run_manifest()` so manifest and trace correlate | High |
| Per-run output | Default to a per-run file under `.metadata/` instead of truncating | Medium |
| Environment attributes | Add `rule.software_environment`, `rule.env_type` to rule spans | Medium |
| `SpanBuilder` | Use it in the emitter or remove it | Low |
| YAML config | Optional `telemetry:` block in benchmark YAML for endpoint config | Low |

## 7. Alternatives Considered

### Alternative 1: Native OpenTelemetry SDK

Use the `opentelemetry-sdk` Python package for span management instead of the custom `TelemetryEmitter`.

- **Pros**: Standard API, less custom code, automatic context propagation.
- **Cons**: Heavyweight dependency, SDK designed for request-scoped instrumentation not batch jobs, harder to emit spans out-of-order (parent before children), less control over OTLP JSON format.
- **Reason not chosen**: The custom emitter gives full control over when spans are emitted (at rule completion, not via context managers) and avoids adding a large transitive dependency to the core package.

### Alternative 2: Snakemake `--report` / native logging

Rely on Snakemake's built-in HTML report and log files for execution records.

- **Pros**: Zero additional code.
- **Cons**: Not machine-readable in a standard format, not compatible with distributed tracing backends, no span hierarchy, log files scattered across `.logs/`.
- **Reason not chosen**: Doesn't address the structured, OTLP-compatible observability requirement.

### Alternative 3: Write logs only (no spans)

Emit structured log records instead of distributed traces.

- **Pros**: Simpler data model.
- **Cons**: Loses the parent-child hierarchy that makes it easy to correlate rule failures to their module and stage context in a trace viewer.
- **Reason not chosen**: The span hierarchy is the key value-add for complex multi-stage benchmarks.

## 8. References

1. [OpenTelemetry Protocol (OTLP) Specification](https://opentelemetry.io/docs/specs/otlp/)
2. [Aspire Dashboard (OTLP backend)](https://learn.microsoft.com/en-us/dotnet/aspire/fundamentals/dashboard/overview)
3. [Jaeger Distributed Tracing](https://www.jaegertracing.io/)
4. [obmon](https://github.com/omnibenchmark/obmon) — consumer of the telemetry stream
5. [Design 007: Output Folder Layout and Runtime Manifest](007-output-layout.md) — `manifest.json` and its `run_id`
6. [PR #333](https://github.com/omnibenchmark/omnibenchmark/pull/333) — implementation
