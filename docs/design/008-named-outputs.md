# 008: Named Outputs and Module/Benchmark Decoupling

[![Status: Draft](https://img.shields.io/badge/Status-Draft-yellow.svg)](https://github.com/omnibenchmark/docs/design)
[![Version: 0.1](https://img.shields.io/badge/Version-0.1-blue.svg)](https://github.com/omnibenchmark/docs/design)

**Authors**: ben
**Date**: 2026-05-05
**Status**: Draft
**Version**: 0.1
**Supersedes**: N/A
**Reviewed-by**: TBD
**Related Issues**: [#325](https://github.com/omnibenchmark/omnibenchmark/issues/325)

## Changes

| Version | Date | Description | Author |
|---------|------|-------------|--------|
| 0.1 | 2026-05-05 | Initial draft | ben |

## 1. Problem Statement

Today the module/benchmark contract is **path-coupled**: the benchmark spec
declares `outputs[].path` (e.g. `"{dataset}.h5ad"`) and modules must produce
files at exactly those paths, reconstructed from `--output_dir` + `--name` +
a filename-suffix naming convention. Renaming a file in the spec requires a
change in every module that produces it.

This breaks module reusability: the same clustering module cannot be used in
two benchmarks that name their outputs differently without forking the module
code.

## 2. Design Goals

- **Decouple module from spec filenames**: modules declare *what* they emit by
  id; the benchmark decides *where on disk* each named output lands.
- **Backward-compatible**: existing v0.4/v0.5 modules keep working unchanged.
- **Handle bags of files**: support outputs that are collections of files
  (fastq libraries, single-cell folder formats) without unknown-cardinality
  Snakemake complexity.
- **Open the door to validation**: design the output contract in a way that
  stage-level validators can be plugged in later.

### Non-Goals

- Removing `--output_dir` / `--name` in this release (deprecated in v0.8).
- Per-file fan-out from directory outputs (checkpoint semantics, deferred).
- Module-side output declaration cross-checks (deferred to module metadata
  reading, stage 3).
- A validation DSL (validators are script-only initially, stage 2).

## 3. Proposed Solution

### 3.1 Output kinds

Extend `IOFile` with an optional `kind` field:

```yaml
outputs:
  - id: rawdata.h5ad
    path: "{dataset}.h5ad"
    kind: file          # default

  - id: rawdata.fastq_collection
    path: "{dataset}.fastq.zip"
    kind: zip           # bag of files, stored uncompressed
```

`kind: zip` uses ZIP64 with no compression (`ZIP_STORED`). Modules are
responsible for writing the zip themselves. Standard libraries in all target
ecosystems support this:

- Python: `zipfile` (stdlib, auto-enables ZIP64)
- R: `zip` CRAN package (not `utils::zip`, which may not support ZIP64)
- Rust: `zip` crate (verify ZIP64 feature)

Zip is preferred over `directory()` for bags of files because it is a single
Snakemake-tracked file, requires no `checkpoint` machinery, and allows
parallel access for downstream readers when entries are stored uncompressed.

### 3.2 Named output contract (API v0.6)

For benchmarks declaring `api_version: 0.6.0`, the runtime passes one
`--output` flag per declared stage output, in addition to the existing
`--output_dir` and `--name` flags (which remain for backward compatibility
through v0.7, dropping in v0.8).

**Multi-output stage:**

```
<entrypoint> \
    --output_dir $OUTPUT_DIR \
    --name <module_id> \
    --output rawdata.h5ad=/abs/path/dataset1.h5ad \
    --output rawdata.clusters_truth=/abs/path/dataset1.clusters_truth.tsv \
    --input_data /abs/path/to/input.h5ad \
    [param flags]
```

**Single-output stage** (drop the `id=` syntax):

```
<entrypoint> \
    --output_dir $OUTPUT_DIR \
    --name <module_id> \
    --output /abs/path/dataset1.h5ad \
    [param flags]
```

Paths passed via `--output` are absolute. Modules that opt into the new
contract write to the exact path given; modules that ignore the new flags and
use `--output_dir` + `--name` continue to work as before.

### 3.3 Snakemake rule output block (v0.6+)

For v0.6 benchmarks the generated Snakemake rule uses Snakemake's own
named-output syntax so ids are tracked by the workflow engine:

```python
rule one_data__datasets__default:
    output:
        rawdata_h5ad="{input}/one-data/datasets/.abc123/{dataset}.h5ad",
        rawdata_clusters_truth="{input}/one-data/datasets/.abc123/{dataset}.clusters_truth.tsv",
```

Output id sanitization (`.` → `_`, etc.) follows the same scheme as the
existing `input_name_mapping`. Collisions after sanitization are a hard
load-time error.

v0.4 and v0.5 rules keep the existing positional `output: "path",` format
unchanged.

### 3.4 Deprecation timeline

| Version | Change |
|---------|--------|
| 0.6 | `--output` flags added (additive). Named-output rule syntax for v0.6 benchmarks. |
| 0.7 | Validator phantom modules (stage 2). |
| 0.8 | `--output_dir` / `--name` dropped. Modules must use `--output` flags. |

## 4. Alternatives Considered

### Alternative 1: `directory()` outputs

Snakemake's `directory(...)` wrapper handles unknown-cardinality outputs.
Rejected because downstream consumers require a `checkpoint` rule and input
functions to fan out over directory contents, adding significant DAG
complexity. ZIP achieves the same result as a single tracked file.

### Alternative 2: Phantom-module-YAML for validators

Declaring validators in benchmark YAML as module-like entries with their own
repository, commit, and env. Rejected in favor of a filesystem convention
(`validators/<stage-id>/` next to `benchmark.yaml`) that lowers to phantom
modules at resolution time — same machinery, zero user-facing YAML overhead.
This is stage 2.

### Alternative 3: DSL for validation rules

A declarative validation mini-language embedded in the `IOFile.validate:`
field. Deferred: too early to know which primitives are useful across real
benchmarks. Stage 2 ships script-only validators; a DSL can emerge from
observed patterns.

## 5. Implementation Plan

### Stage 1 — Named outputs (this PR, API v0.6.0)

1. **`omnibenchmark/model/benchmark.py`**
   - Add `V0_6_0 = "0.6.0"` to `APIVersion`.
   - Add `kind: Literal["file", "zip"] = "file"` to `IOFile`.
   - Delete unused `Module.outputs` field.

2. **`omnibenchmark/model/resolved.py`**
   - Change `ResolvedNode.outputs` from `List[str]` to `Dict[str, str]`
     (id → resolved path template).
   - Add `output_name_mapping: Dict[str, str]` (sanitized key → original id),
     symmetric to the existing `input_name_mapping`.
   - Keep `get_output_list()` returning `list(outputs.values())` for callers
     that only need paths.

3. **`omnibenchmark/benchmark/_graph.py` and `_node.py`**
   - Thread `IOFile.id` through output expansion so ids survive into
     `ResolvedNode.outputs`.

4. **`omnibenchmark/backend/snakemake.py`**
   - `_write_node_rule`: when `api_version >= V0_6_0`, emit named-output
     syntax using sanitized ids; both `kind: file` and `kind: zip` use a
     plain string (zip is a file to Snakemake). v0.4/v0.5 unchanged.
   - `_write_shell`: when `api_version >= V0_6_0`, append `--output` flags
     (single output → positional path; multi → `id=path` per output).
   - Migrate all `node.outputs[N]` integer-index accesses to dict access.
   - Fail hard at codegen if any two output ids sanitize to the same key.

### Stage 2 — Validators (API v0.7.0)

Convention: `validators/<stage-id>/` directory alongside `benchmark.yaml` is
discovered at resolution time and synthesized into a phantom `ResolvedModule`.
The validator entrypoint is invoked per node with `--output id=path` flags
(same contract as the module). Non-zero exit fails the Snakemake DAG.
Generated rule emits a `.validate/<node>.ok` marker; downstream rewires to
depend on the marker.

### Stage 3 — Module metadata + v0.8 decoupling

- Read each module's `omnibenchmark.yaml` at clone time for module-side output
  declarations (`essential: true|false`). Cross-check at benchmark load.
- Drop `--output_dir` / `--name` (v0.8).

### Testing Strategy

- Unit: `tests/backend/test_snakemake_gen.py` — assert v0.6 rule uses
  named-output syntax, single-output emits positional `--output PATH`,
  multi-output emits `--output id=PATH`, `kind: zip` produces plain string
  (no `directory(...)`).
- Back-compat: v0.4/v0.5 fixtures with no `api_version` keep asserting
  positional `output:` format and old shell contract.
- e2e: `tests/e2e/` — small v0.6 module fixture with `kind: file` and
  `kind: zip` outputs.

## 6. References

1. [Issue #325 — flexibility on naming conventions in modules](https://github.com/omnibenchmark/omnibenchmark/issues/325)
2. [ZIP64 extension (PKWARE)](https://pkware.cachefly.net/webdocs/casestudies/APPNOTE.TXT)
3. [Snakemake named outputs](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#named-output-files)
