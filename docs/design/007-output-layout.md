# 007: Output Folder Layout and Runtime Manifest

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
| 0.1 | 2026-02-19 | Initial draft — output layout + runtime manifest | ben |

## 1. Problem Statement

Every `ob run` invocation materialises a tree of files under the configured output directory (default: `out/`). The layout of that tree — which subdirectories exist, what each one contains, and how the files relate to each other — is implicit in the implementation and has never been formally documented.

Additionally, there is no lightweight, structured record of the *execution environment* in which a benchmark was run. The telemetry system (006) captures a rich distributed trace, but requires `--telemetry json` to be enabled. A minimal "who ran this, on what machine, when" record should always be written, regardless of whether telemetry is active.

This document:

1. Formally specifies the output folder layout produced by `ob run`.
2. Defines the `manifest.json` file that records runtime environment metadata for every run.

## 2. Design Goals

- **Documented**: Every directory and file under `out/` has a defined purpose.
- **Stable run identity**: Each invocation of `ob run` is assigned a UUID that appears in the manifest (and in the telemetry trace when telemetry is enabled), so outputs from multiple runs can be correlated.
- **Zero-config provenance**: Host metadata is captured automatically, without any additional flags.
- **Non-intrusive**: No changes to module execution; manifest is written by the orchestrator, not by individual modules.
- **Extensible**: The manifest schema uses a flat JSON object; new fields can be added without breaking readers.

### Non-Goals

- Replacing the telemetry system (006) — the manifest is a complement, not a substitute.
- Capturing per-rule resource usage (CPU time, peak memory) — that belongs in telemetry rule spans.
- Defining remote storage layouts (see 003).

## 3. Output Folder Layout

After a successful (or partial) `ob run`, the output directory has the following structure:

```
out/                               # configurable via --out-dir
│
├── Snakefile                      # generated explicit Snakefile (all wildcards resolved)
│
├── .metadata/                     # provenance and configuration records
│   ├── manifest.json              # run UUID + host metadata (this document, §4)
│   ├── benchmark.yaml             # verbatim copy of the benchmark YAML used for this run
│   └── modules.txt                # list of resolved modules with repo URLs and commits
│
├── .logs/                         # execution logs
│   ├── snakemake_<timestamp>.log  # captured stdout/stderr from the snakemake process
│   └── <rule_name>.log            # per-rule log (stdout/stderr of the module command)
│
├── .modules/                      # checked-out module source trees (git cache worktrees)
│   └── <repo>/<commit>/           # one directory per (repository, commit) pair
│
├── telemetry.jsonl                # OTLP JSON Lines trace (only when --telemetry json)
│
└── <stage_id>/                    # benchmark outputs, one top-level dir per stage
    └── <module_id>/
        └── .<param_hash>/         # deterministic 8-char SHA256 hash of parameters
            └── <output_files>     # as declared in the benchmark YAML outputs[].path
```

### Directory descriptions

| Path | Purpose | Always present |
|------|---------|----------------|
| `Snakefile` | Fully materialised Snakefile; all wildcards resolved to concrete paths. Can be run standalone with `snakemake --cores N`. | Yes |
| `.metadata/` | Human- and machine-readable provenance. Safe to archive with results. | Yes |
| `.metadata/manifest.json` | Run UUID, timestamp, and host metadata (see §4). | Yes |
| `.metadata/benchmark.yaml` | Verbatim copy of the benchmark YAML. | Yes |
| `.metadata/modules.txt` | Human-readable list of modules, their repository URLs, commit hashes, and entrypoints. | Yes |
| `.logs/` | All execution logs. | Yes (created before Snakemake runs) |
| `.logs/snakemake_<ts>.log` | Full stdout/stderr from the Snakemake process. Timestamped to allow multiple runs in the same output dir. | Yes |
| `.logs/<rule>.log` | Per-rule output captured by `tee` inside the shell command. | Only for executed rules |
| `.modules/` | Source trees for resolved modules. Populated during the resolution phase. | Yes |
| `telemetry.jsonl` | OTLP-compatible distributed trace (see 006). | Only with `--telemetry json` |
| `<stage>/` | Stage output directories, nested by `module_id/param_hash/`. | After execution |

### Notes on output nesting

The default nesting strategy (`nested`) mirrors the upstream DAG:

```
<stage_1>/<module_1>/.<hash_1>/<file>
<stage_1>/<module_1>/.<hash_1>/<stage_2>/<module_2>/.<hash_2>/<file>
```

The `<param_hash>` component is prefixed with `.` (dot) so it is hidden by default on Unix systems and clearly distinguished from module IDs.  The hash is the first 8 characters of the SHA256 of the sorted key=value parameter string (see 004 §3.8).

## 4. Runtime Manifest (`manifest.json`)

### 4.1 Purpose

`out/.metadata/manifest.json` is a JSON file written at the start of every `ob run` invocation (before Snakemake executes). It captures:

- A stable **run UUID** that uniquely identifies this invocation.
- Basic **host metadata** (hostname, CPU, RAM, kernel, Python version).
- The **timestamp** at which the run started.

The run UUID correlates with the OTLP `trace_id` in `telemetry.jsonl` when telemetry is enabled, so both files describe the same run.

### 4.2 Schema

```json
{
  "run_id": "550e8400-e29b-41d4-a716-446655440000",
  "ob_version": "0.4.0.post19+g4d4bdce5c.d20260219",
  "snakemake_cmd": ["snakemake", "--snakefile", "Snakefile", "--cores", "8", "--use-conda"],
  "timestamp": "2026-02-19T14:30:00.123456+00:00",
  "hostname": "compute-node-07",
  "platform": "linux",
  "os": "Linux 6.1.0-21-amd64 #1 SMP Debian 6.1.90-1 (2024-05-03)",
  "kernel": "6.1.0-21-amd64",
  "cpu_count": 32,
  "cpu_model": "Intel(R) Xeon(R) Gold 6148 CPU @ 2.40GHz",
  "memory_total_mb": 257024,
  "python_version": "3.11.8",
  "python_executable": "/home/user/.pixi/envs/default/bin/python",
  "gpu_devices": [
    {"index": 0, "name": "NVIDIA RTX 4090", "memory_total_mb": 24576}
  ]
}
```

### 4.3 Field reference

| Field | Type | Description |
|-------|------|-------------|
| `run_id` | string (UUID4) | Unique identifier for this invocation. Matches the OTLP `trace_id` (hex, no dashes) in `telemetry.jsonl` when telemetry is active. |
| `ob_version` | string or null | `omnibenchmark` package version from `importlib.metadata`. Includes git commit and date for development installs (e.g. `0.4.0.post19+g4d4bdce5c.d20260219`). `null` if not installed as a package. |
| `snakemake_cmd` | array of string or null | Exact snakemake argv list as passed to the subprocess (e.g. `["snakemake", "--snakefile", "Snakefile", "--cores", "8"]`). Patched into the manifest by `_run_snakemake` after the command is fully assembled. `null` (field absent) on `--dry` runs. |
| `timestamp` | string (ISO-8601 UTC) | Wall-clock time when the manifest was written (early in `ob run`, before Snakemake starts). |
| `hostname` | string | Hostname of the machine running `ob run`. |
| `platform` | string | Python `sys.platform` value (e.g. `"linux"`, `"darwin"`, `"win32"`). |
| `os` | string | Full OS description from `platform.uname()`: system + release + version. |
| `kernel` | string | Kernel release string (equivalent to `uname -r`). |
| `cpu_count` | integer or null | Logical CPU count visible to the process (`os.cpu_count()`). |
| `cpu_model` | string or null | CPU model name from `/proc/cpuinfo` (Linux) or `sysctl` (macOS). `null` if unavailable. |
| `memory_total_mb` | integer or null | Total physical RAM in MiB from `/proc/meminfo` (Linux) or `sysctl` (macOS). `null` if unavailable. |
| `python_version` | string | Python version string (e.g. `"3.11.8"`). |
| `python_executable` | string | Absolute path to the Python interpreter used to run `ob`. |
| `gpu_devices` | array or null | List of NVIDIA GPUs detected via `nvidia-smi`. Each entry: `{"index": int, "name": str, "memory_total_mb": int}`. `null` if `nvidia-smi` is absent, times out, or returns a non-zero exit code. |

### 4.4 Relationship to telemetry

When `--telemetry json` is active:

- `manifest.json` `run_id` = `xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx` (UUID4 with dashes)
- `telemetry.jsonl` `trace_id` = same UUID in 32-char hex (no dashes)

Conversion: strip dashes from `run_id` to get `trace_id`; insert dashes at positions 8, 12, 16, 20 to convert back.

When telemetry is disabled, `run_id` is a fresh UUID4 generated at the start of `ob run`. It is not present in any other output file, but it uniquely identifies the run in any archival or reporting context.

### 4.5 Stability guarantees

- `run_id` and `timestamp` are written once per `ob run` invocation and never overwritten.
- If `ob run` is re-executed in the same output directory (e.g. to resume after a failure), a new `manifest.json` is written, replacing the previous one. The `run_id` changes for each invocation.
- The schema is additive: new fields may be added in future versions without a version bump. Readers should ignore unknown fields.

## 5. Implementation

### 5.1 `write_run_manifest()` — `omnibenchmark/backend/snakemake_gen.py`

```python
write_run_manifest(output_dir: Path, run_id: Optional[str] = None) -> None
```

- Creates `<output_dir>/.metadata/` if absent.
- Collects system information using `platform`, `sys`, `os`, `/proc/cpuinfo`, `/proc/meminfo`, `sysctl`, and (optionally) `nvidia-smi`.
- Writes `manifest.json` as pretty-printed JSON (2-space indent, trailing newline).
- All system-info collection is wrapped in `try/except`; fields that cannot be read are written as `null`.

### 5.2 Call site — `omnibenchmark/cli/run.py`

`write_run_manifest` is called in `_run_benchmark()`, after `_generate_explicit_snakefile()` completes and the telemetry emitter (if any) has been initialised:

```python
run_id = telemetry_emitter.run_id if telemetry_emitter else None
write_run_manifest(output_dir=out_dir_path, run_id=run_id)
```

### 5.3 `run_id` property — `TelemetryEmitter`

```python
@property
def run_id(self) -> str:
    """Return the run UUID (same as OTLP trace_id, formatted as standard UUID string)."""
    t = self._trace_id  # 32-char hex
    return f"{t[0:8]}-{t[8:12]}-{t[12:16]}-{t[16:20]}-{t[20:32]}"
```

## 6. Alternatives Considered

### Alternative 1: Write manifest only when telemetry is active

Keep `manifest.json` as a telemetry-only artifact.

- **Pros**: Simpler — one less always-on side effect.
- **Cons**: Loses provenance for the majority of runs (telemetry is opt-in). A run without telemetry produces no machine-readable record of when or where it ran.
- **Reason not chosen**: The manifest is lightweight and provides value independent of telemetry.

### Alternative 2: Embed host metadata in `modules.txt`

Append machine info to the existing `modules.txt` file.

- **Pros**: No new file.
- **Cons**: `modules.txt` is a human-readable text file, not machine-readable. Mixing host metadata with module provenance makes both harder to parse. The run UUID has no natural home in `modules.txt`.
- **Reason not chosen**: A dedicated JSON file is cleaner and more extensible.

### Alternative 3: Use `psutil` for hardware info

Use the `psutil` library for portable CPU/RAM collection instead of reading `/proc` and `sysctl` directly.

- **Pros**: Cross-platform, richer metrics.
- **Cons**: Adds a dependency. The current approach covers Linux and macOS (the two primary platforms) with zero new dependencies, gracefully falling back to `null` on others.
- **Reason not chosen**: Avoiding new dependencies is preferred for core functionality. Can be revisited if Windows support becomes a requirement.

## 7. References

1. [Design 004: YAML Specification](004-yaml-specification.md) — parameter hash algorithm
2. [Design 006: Telemetry](006-telemetry.md) — OTLP trace_id format and span hierarchy
3. [OpenTelemetry Trace ID spec](https://opentelemetry.io/docs/specs/otel/overview/#trace-id)
