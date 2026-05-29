# 004: Omnibenchmark YAML Specification (v0.1)

[![Status: Under Review](https://img.shields.io/badge/Status-UnderReview-yellow.svg)](https://github.com/omnibenchmark/docs/design)
[![Version: 0.1](https://img.shields.io/badge/Version-0.3-blue.svg)](https://github.com/omnibenchmark/docs/design)

**Authors**: ben
**Date**: 2025-01-20
**Status**: Under Review
**Version**: 0.1
**Supersedes**: N/A
**Reviewed-by**: daninci
**Related Issues**: #283

## Changes

| Version | Date | Description | Author |
|---------|------|-------------|--------|
| 1       | 2026-01-20 | Initial specification for v0.4 | ben |
| 2       | 2026-03-31 | Add resource allocation (Section 8) | ben |
| 3       | 2026-05-06 | Add provenance metadata (Section 9): canonical_url, derived_from, subset_of | ben |

## 1. Problem Statement

Omnibenchmark requires a formal specification for its YAML-based benchmark definition format. This specification serves as the canonical reference for:

- Benchmark authors defining computational pipelines
- Tool implementers building parsers and compilers
- Backend generators targeting workflow engines (Snakemake, etc.)

Without a formal specification, the format's semantics are implicit in the implementation, leading to ambiguity and maintenance challenges.

## 2. Design Goals

- **Declarative**: Separate benchmark logic from execution engine details
- **Modular**: Enable reusable components with clear interfaces
- **Reproducible**: Explicit versioning and environment specifications
- **Scalable**: Automatic cartesian product expansion across modules
- **Backend-agnostic**: Compile to multiple workflow engines

### Non-Goals

- Runtime execution details (handled by backends)
- Advanced features like conditional execution, parameterized modules, specialized data modules, or dynamic workflows (future versions)
- Non-YAML formats (JSON Schema, etc.)

## 3. Specification

### 3.1 Document Structure

A valid omnibenchmark YAML document has the following top-level structure:

```yaml
id: <string>                           # Required: unique identifier
benchmarker: <string>                  # Required: author/organization
version: <string>                      # Required: semantic version
description: <string>                  # Optional: human-readable description
api_version: <string>                  # Optional: spec version (default: "0.4")
software_backend: <string>             # Optional: default backend (e.g., "conda")
software_environments: <object>        # Optional: named environments
stages: <array>                        # Required: pipeline stages (≥1)
provenance:                            # Optional: lineage metadata
  canonical_url: <string>             #   URL of the authoritative publication
  derived_from: <string>              #   identifier/URL of the upstream benchmark
  subset_of: <string>                 #   summary hash of the parent benchmark
```

#### Required Fields

- `id`: Alphanumeric identifier (underscores allowed)
- `benchmarker`: Author or organization name
- `version`: Semantic version string
- `stages`: Array of stage definitions (minimum 1)

#### Optional Fields

- `description`: Human-readable benchmark description
- `api_version`: Specification version (default: `"0.3"`)
- `software_backend`: Default execution backend
- `software_environments`: Named environment definitions
- `provenance`: Lineage metadata (see Section 9)

### 3.2 Software Environments

Software environments define execution contexts for modules.

#### Host Environment

Direct execution on host system without isolation:

```yaml
software_environments:
  host:
    description: "Direct host execution"
```

This environment is useful for direct execution on the host system without any isolation, and as such it should be used with caution and discouraged for production benchmarks.

#### Container Environment

Execution within containers (Docker, Apptainer, etc.):

```yaml
software_environments:
  my_container:
    description: "Containerized environment"
    apptainer: "docker://ubuntu:22.04"
```

#### Conda Environment

Execution within conda environments:

```yaml
software_environments:
  my_conda:
    description: "Conda environment"
    conda: "environment.yml"  # Path to conda environment file
```

#### Environment Resolution

When a module references `software_environment: "env_name"`:

1. Lookup `env_name` in `software_environments`
2. Validation error if not found
3. Backend receives resolved environment for directive generation

### 3.3 Stages

Stages represent logical phases in the pipeline.

```yaml
stages:
  - id: <string>                       # Required: unique stage identifier
    modules: <array>                   # Required: module list (≥1)
    inputs: <array>                    # Optional: input dependencies
    outputs: <array>                   # Optional: output declarations
    provides: <array>                  # Optional: lineage labels (api 0.6+)
```

`provides:` is a list of label names this stage advertises to downstream
modules; downstream modules can gate on these labels via `requires:`. See
[008-filtering.md §3.5](./008-filtering.md) for details.

#### Stage Ordering

- Stages are processed in document order
- Dependencies resolved via input/output declarations
- Stage N consumes outputs from Stage N-1

### 3.4 Modules

Modules are concrete implementations within stages.

```yaml
modules:
  - id: <string>                       # Required: unique module ID
    software_environment: <string>     # Required: environment reference
    repository: <object>               # Required: source location
    name: <string>                     # Optional: human-readable name
    parameters: <array>                # Optional: module parameters
    exclude: <array>                   # Optional: exclusion list
    requires_capabilities: <array>     # Optional: required host capabilities
```

#### Required Fields

- `id`: Unique identifier within stage (used as wildcard value)
- `software_environment`: Reference to defined environment
- `repository`: Source code specification

#### Capability Gating

`requires_capabilities` declares host capabilities the module needs in order
to execute (e.g. `gpu`, `large_mem`). At run time, modules whose required set
is not a subset of the capabilities provided via `ob run --capability <name>`
(repeatable) are silently pruned from the resolved DAG. Their outputs are
never registered, so downstream stages naturally lose those branches without
any explicit exclusion logic.

```yaml
modules:
  - id: pca_gpu
    software_environment: cuda
    repository:
      url: https://github.com/example/pca-gpu
      commit: abc1234
    requires_capabilities: [gpu]
```

Capabilities are deliberately scoped to *host* requirements. Free-form tags
("experimental", "phase2", …) are out of scope; if a value begs for
free-form text it should be a parameter, a stage split, or a separate
benchmark. See [008-filtering.md](./008-filtering.md) for the broader
filtering design.

#### Repository Specification

```yaml
repository:
  url: <string>      # Repository URL (git, local path, etc.)
  commit: <string>   # Commit hash, tag, or branch
```

### 3.5 Parameters

Parameters are passed to module execution as GNU-style command-line arguments.

```yaml
parameters:
  - key1: value1
    key2: value2
  - key3: value3
```

Resolves to: `--key1 value1 --key2 value2 --key3 value3`

#### Reserved Parameters

The system automatically provides:

- `--name`: Module ID for current execution (e.g., `--name D1`)
- `--output_dir`: Output directory path (e.g., `--output_dir data/D1/29b6dbbe`)

Modules MUST NOT declare these in their YAML configuration.

#### Parameter Values

- **String**: Passed directly (`evaluate: "1+1"` → `--evaluate 1+1`)

#### Wildcard-Based Resolution

Parameters are mapped to wildcard values:

```yaml
modules:
  - id: D1
    parameters:
      - evaluate: "1+1"
  - id: D2
    parameters:
      - evaluate: "2+2"
```

Backend mapping (Snakemake example):

```python
params:
    evaluate=lambda wildcards: {"D1": "1+1", "D2": "2+2"}[wildcards.dataset]
```

### 3.6 Inputs and Outputs

#### Output Declarations

Stages declare what they produce:

```yaml
outputs:
  - id: data.raw                       # Output identifier
    path: "{dataset}_data.json"        # Path template with wildcards
```

- `id`: Used for referencing in downstream stages
- `path`: Template with `{wildcard}` placeholders. Here, {dataset} is a special case, which is replaced with the module ID of the initial stage.

#### Input Declarations

Stages declare dependencies on upstream outputs:

```yaml
inputs:
  - data.raw                           # Simple reference
```

#### Input Resolution

The compiler maintains an output registry:

1. Register outputs: `data.raw → "data/{dataset}/{dataset}_data.json"`
2. Resolve inputs by lookup
3. Pass to modules: `--data.raw data/{dataset}/{dataset}_data.json`

### 3.7 Path Strategy

#### Default Nesting

Hierarchical structure:

```
<stage_id>/<wildcard_pattern>/<params_hash>/<output_file>
```

Example:
- Stage: `data`, Module: `D1`, Output: `{dataset}_data.json`
- Resolved: `data/{dataset}/29b6dbbe/{dataset}_data.json`

#### Parameter Hash

A deterministic hash distinguishes parameter combinations:

**Algorithm:**
1. Extract all parameters from module
2. Sort key-value pairs lexicographically by key
3. Serialize to canonical form: `"key1=value1,key2=value2,..."`
4. Apply SHA256 hash
5. Truncate to first 8 characters

**Example:**

```python
# Parameters: {evaluate: "1+1", size: "100"}
# Sorted: [("evaluate", "1+1"), ("size", "100")]
# Canonical: "evaluate=1+1,size=100"
# SHA256: "f4a3b2c1d5e6f7a8b9c0d1e2f3a4b5c6..."
# Result: "f4a3b2c1"
```

**Special cases:**
- Empty parameters → hash of empty string
- List values → comma-joined before hashing

### 3.8 Wildcard System

#### Wildcard Naming

Wildcard names are extracted from output path templates:

- `{dataset}_data.json` → wildcard name: `dataset`
- `result_{method}.txt` → wildcard name: `method`

**Fallback:** If no wildcard in template, use lowercased stage ID.

**Conflict resolution:** If extracted wildcard already inherited, use stage ID.

#### Wildcard Values

Module IDs become wildcard values:

```yaml
stages:
  - id: data
    modules:
      - id: D1
      - id: D2
    outputs:
      - path: "{dataset}_data.json"
```

Results in: `{dataset: [D1, D2]}`

#### Wildcard Propagation

Wildcards accumulate through pipeline, creating cartesian products:

```
Stage 1 (data):    {dataset: [D1, D2]}
Stage 2 (methods): {dataset: [D1, D2], method: [M1, M2]}
Stage 3 (metrics): {dataset: [D1, D2], method: [M1, M2], metric: [R1, R2]}
```

Downstream stages inherit wildcards from upstream dependencies.

#### Cartesian Exclusions

Modules can exclude specific combinations:

```yaml
modules:
  - id: D2
    exclude: [M2]  # Exclude combination dataset=D2 AND method=M2
```

## 4. Complete Example

### Simple Multi-Stage Pipeline

```yaml
id: MethodBenchmark
version: "1.0"
benchmarker: "Benchmark Team"
api_version: "0.3"

software_environments:
  host:
    description: "Host execution"

stages:
  # Stage 1: Data generation
  - id: data
    modules:
      - id: D1
        software_environment: "host"
        repository:
          url: "bundles/data.bundle"
          commit: "abc123"
        parameters:
          - n: "100"
      - id: D2
        software_environment: "host"
        repository:
          url: "bundles/data.bundle"
          commit: "abc123"
        parameters:
          - n: "1000"
    outputs:
      - id: data.raw
        path: "{dataset}_data.json"

  # Stage 2: Method execution
  - id: methods
    modules:
      - id: M1
        software_environment: "host"
        repository:
          url: "bundles/method.bundle"
          commit: "def456"
        parameters:
          - algo: "fast"
      - id: M2
        software_environment: "host"
        repository:
          url: "bundles/method.bundle"
          commit: "def456"
        parameters:
          - algo: "accurate"
    inputs:
      - data.raw
    outputs:
      - id: methods.result
        path: "{dataset}_{method}_result.json"
```

### Compiled Snakemake (conceptual)

```python
# Cartesian product: D1×M1, D1×M2, D2×M1, D2×M2
rule all:
    input:
        expand("methods/{dataset}/{method}/29b6dbbe/{dataset}_{method}_result.json",
               dataset=["D1", "D2"],
               method=["M1", "M2"])

rule data:
    output:
        "data/{dataset}/a1b2c3d4/{dataset}_data.json"
    params:
        n=lambda wildcards: {"D1": "100", "D2": "1000"}[wildcards.dataset]
    shell:
        "data_gen --name {wildcards.dataset} --output_dir data/{wildcards.dataset}/a1b2c3d4 --n {params.n}"

rule methods:
    input:
        "data/{dataset}/a1b2c3d4/{dataset}_data.json"
    output:
        "methods/{dataset}/{method}/29b6dbbe/{dataset}_{method}_result.json"
    params:
        algo=lambda wildcards: {"M1": "fast", "M2": "accurate"}[wildcards.method]
    shell:
        "method --name {wildcards.method} --output_dir methods/{wildcards.dataset}/{wildcards.method}/29b6dbbe --data.raw {input} --algo {params.algo}"
```

## 5. Implementation Notes

### Current Limitations

- **Nested parameters**: Object-valued parameters not fully supported
- **Cartesian exclusions**: Simplified implementation needs enhancement
- **Wildcard extraction**: Limited conflict resolution logic

### Backend Requirements

Backends must support:

1. **Cartesian product expansion**: Generate all wildcard combinations
2. **Parameter mapping**: Map module parameters to wildcard values
3. **Path templating**: Substitute wildcards in path templates
4. **Environment directives**: Translate environment specs to backend syntax
5. **Dependency resolution**: Wire inputs to upstream outputs

## 7. Resource Allocation

> **Status:** Implemented.

### 7.1 Overview

Stages, modules, and metric collectors can declare resource requirements via an optional `resources:` block. These are translated directly into Snakemake `resources:` directives, enabling efficient parallel scheduling.

### 7.2 Schema

```yaml
resources:
  cores: <integer>      # Optional: CPU cores per task
  mem_mb: <integer>     # Optional: memory in MiB
  disk_mb: <integer>    # Optional: disk space in MiB
  runtime: <integer>    # Optional: expected runtime in minutes
```

All fields are optional, but at least one must be present if a `resources:` block is declared.

| Field | Type | Unit | Description |
|-------|------|------|-------------|
| `cores` | integer | — | Logical CPU cores to allocate per task |
| `mem_mb` | integer | MiB | Physical memory to allocate per task |
| `disk_mb` | integer | MiB | Disk space required per task |
| `runtime` | integer | minutes | Expected wall-clock runtime (used by Snakemake to prioritise long tasks) |

### 7.3 Placement

`resources:` can appear at three levels:

**Stage level** — applies to all modules in the stage:

```yaml
stages:
  - id: clustering
    resources:
      cores: 20
      mem_mb: 16000
    modules:
      - id: fastcluster
      - id: sklearn
```

**Module level** — overrides the stage default for that module:

```yaml
stages:
  - id: clustering
    resources:
      cores: 20
      mem_mb: 16000
    modules:
      - id: fastcluster
        # inherits stage resources
      - id: lightweight
        resources:
          cores: 4
          mem_mb: 4000
```

**Metric collector level**:

```yaml
metric_collectors:
  - id: plotting
    resources:
      cores: 4
      mem_mb: 8000
```

### 7.4 Resolution Hierarchy

Resources are resolved with **most-specific-wins** precedence:

1. Module-level `resources:` (highest priority)
2. Stage-level `resources:`
3. Global default: `cores=2` (applied when no `resources:` block is present anywhere in the chain)

### 7.5 Snakemake Mapping

Each non-`null` resource field is emitted as a Snakemake `resources:` entry. The `cores` field additionally sets the `threads:` directive so Snakemake accounts for CPU consumption in its scheduler.

```python
rule clustering_fastcluster_abc12345:
    threads: 20
    resources:
        cores=20,
        mem_mb=16000,
    shell: ...
```

Fields absent from the `resources:` block are not emitted (Snakemake uses its own defaults for unspecified resources).

### 7.6 Controlling parallelism at runtime

Pass `--resources` to the generated Snakemake invocation to cap total consumption:

```bash
# Cap to 100 cores and 128 GiB RAM across all concurrent rules
snakemake --cores 100 --resources mem_mb=131072
```

Snakemake's greedy scheduler will run as many rules in parallel as fit within the declared resource pool.

### 7.7 Complete example

```yaml
id: ResourceExample
version: "1.0"
benchmarker: "Benchmark Team"

software_environments:
  host:
    description: "Host execution"

stages:
  - id: data
    # No resources: global default applies (cores=2)
    modules:
      - id: generator
        software_environment: host
        repository: { url: "...", commit: "abc" }

  - id: clustering
    resources:
      cores: 20
      mem_mb: 16000
    modules:
      - id: fastcluster
        software_environment: host
        repository: { url: "...", commit: "def" }

      - id: lightweight
        resources:
          cores: 4
          mem_mb: 4000
        software_environment: host
        repository: { url: "...", commit: "def" }

metric_collectors:
  - id: plotting
    resources:
      cores: 4
      mem_mb: 8000
    software_environment: host
    repository: { url: "...", commit: "ghi" }
    inputs: [clustering.result]
    outputs:
      - id: plot
        path: "plot.html"
```

## 9. Provenance Metadata

The optional `provenance` block records lineage relationships between benchmarks. All fields are informational — the runtime does not act on them today, but they are preserved in the parsed model for tooling and future use.

```yaml
provenance:
  canonical_url: "https://omnibenchmark.org/b/clustering-benchmark"
  derived_from: "clustering-benchmark"
  subset_of: "a1b2c3d4"
```

### Fields

| Field | Type | Description |
|-------|------|-------------|
| `canonical_url` | string (URL) | Where the authoritative, published version of this benchmark lives. Useful when a benchmark is mirrored or forked across organisations. |
| `derived_from` | string | Identifier or URL of the upstream benchmark this one was forked or adapted from. Intended to be populated automatically by a future `ob fork` command. |
| `subset_of` | string | Summary hash of the parent benchmark of which this is a strict subset (e.g. fewer stages, modules, or parameter grid). The hash should match the value returned by `Benchmark.summary_hash()` once that method is public. |

### Rationale

These fields answer three common questions:

- **Where is the canonical version?** (`canonical_url`) — multiple groups may publish mirrors or slight variants of the same benchmark. This field points back to the source of truth.
- **Where did this come from?** (`derived_from`) — tracks intellectual lineage. When/if a `fork` command is added, it will write this field automatically.
- **Is this a strict subset?** (`subset_of`) — a reduced benchmark (e.g. a quick-run variant with fewer methods) can declare which full benchmark it was derived from, enabling tools to aggregate or compare results across the two. It also allows for automated filtering and grouping of benchmarks by their subset relationship. For instance, a web editor can offer a wizard to select subsets of a benchmark for "slicing and dicing" the plan themselves, and export the subset in a standard compact notation.

### Notes on `subset_of` hashing

`subset_of` should be the value returned by `Benchmark.summary_hash()` on the parent benchmark. This is the same hash used to tag storage artifacts (see design/003-storage.md, Section 3 — "HASH of benchmark yaml"). It covers the execution-relevant fields only: `id`, `software_backend`, `software_environments`, `stages`, and `metric_collectors`. Descriptive fields (`benchmarker`, `authors`, `description`), the `version` string (already encoded separately in the artifact version tag as `VERSION-HASH`), storage configuration, `api_version`, and `provenance` itself are all excluded so that purely administrative edits do not invalidate artifact lineage.

The hash follows the same canonicalize-then-SHA256 convention as parameter set hashing: the execution fields are serialised to canonical JSON (`json.dumps(..., sort_keys=True)`) and the full 64-character hex digest is the output. The short 8-character form (first 8 hex chars) matches the `HASH` component of the `VERSION-HASH` artifact tag described in 003-storage.md.

## 8. References

1. [Omnibenchmark Documentation](https://docs.omnibenchmark.org)
2. [Snakemake Documentation](https://snakemake.readthedocs.io)
3. [Semantic Versioning](https://semver.org)
