# 006: Benchmark Execution Telemetry

[![Status: Draft](https://img.shields.io/badge/Status-Draft-yellow.svg)](https://github.com/omnibenchmark/docs/design)
[![Version: 0.1](https://img.shields.io/badge/Version-0.1-blue.svg)](https://github.com/omnibenchmark/docs/design)

**Authors**: ben
**Date**: 2026-02-19
**Status**: Draft
**Version**: 0.1
**Supersedes**: N/A
**Reviewed-by**: TBD
**Related Issues**: TBD

## Changes

| Version | Date | Description | Author |
|---------|------|-------------|--------|
| 0.1 | 2026-02-19 | Initial draft | ben |

## 1. Problem Statement

Benchmark runs are long-lived, distributed computations. Understanding what happened — which rules ran, how long they took, what failed and why — currently requires reading Snakemake logs, per-rule log files, and stdout. There is no structured, machine-readable record of a benchmark execution that downstream tools (dashboards, CI systems, result databases) can consume.

Additionally, the current `ob run` execution model produces no progress signal compatible with distributed observability infrastructure. Reproducibility and auditability of benchmark runs require richer provenance than what Snakemake natively provides.

### Current State

A telemetry subsystem is already partially implemented (`omnibenchmark/telemetry/`):

- **What works:** OTLP-compatible JSON Lines output, a four-level span hierarchy (benchmark → stage → module → rule), per-rule stdout/stderr capture, relay script to forward to Aspire/Jaeger.
- **What's incomplete:** No tests, an unused `SpanBuilder` class, no automatic output to a standard location, no YAML-level configuration, telemetry is coupled to the explicit Snakefile execution path.
- **Key dependency:** Telemetry requires the *explicit Snakefile* (`_generate_explicit_snakefile()` in `run.py`) — the fully materialised DAG where all wildcards are resolved to concrete paths before Snakemake runs. It cannot work with Snakemake's native wildcard expansion because the span hierarchy must be known before execution starts.

## 2. Design Goals

- **Structured execution record**: Every benchmark run produces a machine-readable trace — a hierarchical record of what ran, when, how long, and whether it succeeded.
- **OTLP compatibility**: Output conforms to the OpenTelemetry Protocol (OTLP) so it can be forwarded to any standard tracing backend (Jaeger, Tempo, Aspire Dashboard, etc.) without custom tooling.
- **Zero-config baseline**: Telemetry works without any external infrastructure. The default output is a local NDJSON file alongside the benchmark outputs.
- **Non-intrusive**: Telemetry is opt-in at the CLI level. Benchmarks that don't use `--telemetry` are unaffected.
- **Explicit-Snakefile dependency is explicit**: The design acknowledges and documents that telemetry is only available when using the explicit Snakefile execution path.

### Non-Goals

- Replacing Snakemake's native logging or `--report`.
- Real-time streaming dashboards (possible via relay, but not a built-in requirement).
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

### 3.2 Span Hierarchy

A benchmark run is represented as a four-level span tree:

```
benchmark_run                        # root: one per ob run invocation
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

**Root span** attributes:
- `benchmark.name`, `benchmark.version`, `benchmark.author`
- `benchmark.total_rules`, `benchmark.software_backend`, `benchmark.cores`

**Rule span** attributes:
- `rule.name`, `rule.stage_id`, `rule.module_id`, `rule.node_id`, `rule.param_id`
- `rule.inputs` (list of paths), `rule.outputs` (list of paths)
- `rule.parameters` (JSON-encoded parameter dict)

All spans carry a common `trace_id` so the full run can be correlated across the hierarchy.

### 3.3 CLI Interface

```
ob run benchmark.yaml [--telemetry json] [--telemetry-output PATH]
```

- `--telemetry json`: Enable telemetry output. Writes to `<out_dir>/telemetry.jsonl` by default. Rich progress UI remains active.
- `--telemetry-output PATH`: Override the default output location.

### 3.4 Relay to External Backends

A relay script (`scripts/telemetry-relay.py`) converts NDJSON to OTLP protobuf and forwards to a gRPC endpoint:

```bash
# Stream live telemetry to Aspire Dashboard (telemetry written to out/telemetry.jsonl by default)
ob run benchmark.yaml --telemetry json &
tail -f out/telemetry.jsonl | python scripts/telemetry-relay.py --endpoint localhost:4317

# Or via convenience script (uses default path)
ob run benchmark.yaml --telemetry json &
./stream_telemetry.sh
```

No OpenTelemetry SDK dependency is required by the core omnibenchmark package. The relay script is optional and listed as a dev dependency.

### 3.5 Dependency on Explicit Snakefile

Telemetry requires the explicit Snakefile execution path. This is the only supported path for `ob run` today, but the dependency is worth documenting:

- **Why:** Telemetry initialises the full span hierarchy (benchmark → stage → module → rule) from the resolved node list *before* Snakemake executes. This is only possible when the complete DAG is materialised upfront (i.e. wildcards resolved, all node IDs known).
- **Consequence:** If a future backend uses native Snakemake wildcard expansion (not explicit rules), telemetry as specified here cannot work without a different integration approach.
- **Constraint captured:** `ob run` always uses `_generate_explicit_snakefile()`. Any backend that bypasses explicit rule generation loses telemetry support.

## 4. Current Implementation Status

| Component | Status | Notes |
|-----------|--------|-------|
| `TelemetryEmitter` | Working | Core emitter, ~740 lines, no external deps |
| OTLP JSON output | Working | Spec-compliant NDJSON |
| `--telemetry` / `--telemetry-output` CLI flags | Working | Integrated in `ob run` |
| Span hierarchy (benchmark/stage/module/rule) | Working | Four levels, full trace_id correlation |
| Per-rule stdout/stderr capture | Working | Via shell `tee` + log file |
| `scripts/telemetry-relay.py` | Working | NDJSON → OTLP gRPC |
| `SpanBuilder` class | Removed | Dead code deleted |
| `use_telemetry_stdout` branch | Removed | Dead code deleted |
| Automatic output location | Done | Defaults to `<out_dir>/telemetry.jsonl` |
| Tests | Missing | No unit or integration tests |
| YAML configuration | Missing | No benchmark-level telemetry config |
| Software environment in spans | Missing | Environment details absent from rule spans |

## 5. Known Issues and Gaps

**No tests.** The telemetry module has no unit or integration tests. Regression risk is high given the size (~740 lines) of `TelemetryEmitter`.

**Environment metadata missing.** Rule spans don't include which software environment was used, whether it was a container or conda env, or what the resolved image/env path was.

## 6. Planned Improvements

| Item | Description | Priority |
|------|-------------|----------|
| Unit tests | Test span emission, hierarchy construction, OTLP format | High |
| Environment attributes | Add `rule.software_environment`, `rule.env_type` to rule spans | Medium |
| Run ID | Stable, user-visible run identifier propagated through spans | Medium |
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
