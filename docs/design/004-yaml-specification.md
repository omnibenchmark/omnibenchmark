# 004: Omnibenchmark YAML Specification (v0.5)

[![Status: Under Review](https://img.shields.io/badge/Status-Accepted-green.svg)](https://github.com/omnibenchmark/docs/design)
[![Version: 0.4](https://img.shields.io/badge/Version-0.3-blue.svg)](https://github.com/omnibenchmark/docs/design)

**Authors**: ben,
**Date**: 2025-01-20
**Status**: Under Review
**Version**: 0.6
**Supersedes**: N/A
**Reviewed-by**: TBD
**Related Issues**: TBD

## Changes

| Version | Date | Description | Author |
|---------|------|-------------|--------|
| 0.1 | 2026-01-20 | Initial specification for v0.4.0 | ben |
| 0.2 | 2026-02-09 | Add named entrypoints (Section 3.5) | ben |
| 0.3 | 2026-02-19 | Add gather stages and path templates (Section 7) | ben |
| 0.4 | 2026-02-19 | Add resource allocation (Section 8) | ben |

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
- **Backend-agnostic**: Compile to multiple workflow engines in the future, with primary focus on Snakemake for the time being

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
```

#### Required Fields

- `id`: Alphanumeric identifier (underscores allowed)
- `benchmarker`: Author or organization name
- `version`: Semantic version string
- `stages`: Array of stage definitions (minimum 1)

#### Optional Fields

- `description`: Human-readable benchmark description
- `api_version`: Specification version, matching omnibenchmark version (default to the version that produced the config, e.g. `"0.4"`)
- `software_backend`: Default execution backend
- `software_environments`: Named environment definitions

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
```

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
```

#### Required Fields

- `id`: Unique identifier within stage (used as wildcard value)
- `software_environment`: Reference to defined environment
- `repository`: Source code specification

#### Repository Specification

```yaml
repository:
  url: <string>            # Required: Repository URL (git, local path, etc.)
  commit: <string>         # Required: Commit hash, tag, or branch
  entrypoint: <string>     # Optional: named entrypoint key (default: "default")
```

### 3.5 Entrypoints

Each module repository contains an `omnibenchmark.yaml` that declares one or more named entrypoints:

```yaml
# Module's omnibenchmark.yaml
entrypoints:
  default: "run.py"
  preprocess: "preprocess.py"
  validate: "scripts/validate.R"
```

#### Default Behavior

When a module's `repository` block in the benchmark plan omits the `entrypoint` field, the system uses the `default` key from the module's `omnibenchmark.yaml`. This preserves full backward compatibility.

#### Named Entrypoints

A benchmark plan can select a specific entrypoint by name. This allows a single repository to expose multiple scripts without requiring separate modules or multiple dispatch:

```yaml
modules:
  - id: preprocess_step
    software_environment: "host"
    repository:
      url: "bundles/pipeline.bundle"
      commit: "abc123"
      entrypoint: "preprocess"       # uses entrypoints.preprocess

  - id: validate_step
    software_environment: "host"
    repository:
      url: "bundles/pipeline.bundle"
      commit: "abc123"
      entrypoint: "validate"         # uses entrypoints.validate
```

#### Resolution Rules

1. If `entrypoint` is omitted or not specified, use `"default"`
2. Look up the key in the module's `omnibenchmark.yaml` `entrypoints` dict
3. If the key is not found, resolution fails with an error listing available entrypoints
4. An empty or whitespace-only `entrypoint` value is a validation error
5. Legacy `config.cfg` modules only support the `"default"` entrypoint; requesting a named entrypoint from a `config.cfg` module is an error

### 3.6 Parameters

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
- **List**: Joined with commas (`items: [a, b, c]` → `--items a,b,c`)

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

### 3.7 Inputs and Outputs

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

### 3.8 Path Strategy

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

### 3.9 Wildcard System

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

## 6. References

1. [Omnibenchmark Documentation](https://docs.omnibenchmark.org)
2. [Snakemake Documentation](https://snakemake.readthedocs.io)
3. [Semantic Versioning](https://semver.org)

---

## 7. Gather Stages (v0.3+)

> **Status:** Implemented (phases 1–6 complete as of 2026-02-08).
> This section documents the `provides` / `gather` extension to the stage model.

### 7.1 Motivation

The pipeline has two execution patterns:

- **Map stages** (`stages`): each module runs once per input combination (cartesian product). First-class, composable.
- **Reduce/gather stages** (`metric_collectors`): a module runs once and receives *all* outputs from referenced stages. Second-class, not composable, not chainable.

Gather stages generalise the reduce pattern: any stage can gather all outputs from any set of provider stages, and gather stage outputs participate in the DAG like regular stage outputs.

### 7.2 New YAML Keywords

Two new keywords extend the stage model:

- **`provides: {label: output_id, ...}`** on a stage — a map from label name to the output ID that satisfies it. Downstream stages can `gather` by label, and downstream map stages get `{label}` as a path template variable resolving to the provider's module ID.
- **`gather: label`** as an entry in a stage's `inputs` list — collects the declared output (by output_id) from every resolved node of all stages that `provides` that label.

```yaml
stages:
  - id: data
    provides:
      dataset: data.raw              # label "dataset" is satisfied by output "data.raw"
    modules:
      - id: D1
        ...
      - id: D2
        ...
    outputs:
      - id: data.raw
        path: "{dataset}_data.json"     # {dataset} = own module ID

  - id: methods_fast
    provides:
      method: methods_fast.result    # label "method" is satisfied by "methods_fast.result"
    inputs: [data.raw]
    modules:
      - id: M1
        ...
    outputs:
      - id: methods_fast.result
        path: "{dataset}_{method}_result.json"

  - id: methods_accurate
    provides:
      method: methods_accurate.result
    inputs: [data.raw]
    modules:
      - id: M2
        ...
    outputs:
      - id: methods_accurate.result
        path: "{dataset}_{method}_result.json"

  - id: summary
    modules:
      - id: S1
        ...
    inputs:
      - gather: method               # collects from methods_fast + methods_accurate
    outputs:
      - id: summary.report
        path: "report.json"
```

### 7.3 Semantics of `provides`

`provides: {label: output_id}` on a stage has two roles:

**Role 1: Gather contract.** Any downstream stage with `gather: label` in its inputs collects the specified output path (by `output_id`) from every resolved node of all stages providing that label. Only the named output is gathered — if a stage has multiple outputs, only the one declared in `provides` is exposed to gatherers.

**Role 2: Template variable.** Any downstream map stage that transitively depends on a providing stage gets `{label}` as a template variable in its output paths. The variable resolves to the module ID of the ancestor that provided that label. For the providing stage itself, `{label}` resolves to its own module ID (self-reference).

`provides` labels live in a separate namespace from output IDs. A stage can have both `provides: {dataset: data.raw}` and `outputs: [{id: dataset.raw, ...}]` without conflict.

Multiple stages can provide the same label. Map stages still reference specific output IDs via `inputs:`; the `provides` label is only used for `gather:` and template variable propagation.

### 7.4 Semantics of `gather:`

`gather: label` in a stage's inputs means:

1. Find every stage that `provides: {label: output_id}`
2. Collect that `output_id`'s concrete output path from every resolved node of those stages
3. Pass all of them as inputs to the gather stage's single rule invocation
4. The gather stage runs **once** per (module × params) combination in the gather stage itself — not once per upstream combination

**Constraints:**
- A stage with `gather:` inputs cannot also have regular `inputs:` (either map or gather, not both)
- A gather stage's output path must not contain provides-derived template variables (it produces a single file per module × params combination)
- A gather stage can have multiple modules and parameters; each combination produces one gather rule that receives all upstream outputs

### 7.5 Composability

Gather stage outputs are regular stage outputs and participate in the DAG:

```yaml
# Gather → Map
- id: summary
  inputs:
    - gather: method
  outputs:
    - id: summary.report
      path: "report.json"

- id: postprocess
  inputs: [summary.report]             # single file, maps once
  outputs:
    - id: postprocess.result
      path: "postprocessed.json"
```

```yaml
# Gather → Gather (chained)
- id: first_gather
  provides:
    intermediate: first_gather.output
  inputs:
    - gather: method
  outputs:
    - id: first_gather.output
      path: "intermediate.json"

- id: second_gather
  inputs:
    - gather: intermediate
  outputs:
    - id: final.output
      path: "final.json"
```

### 7.6 Snakemake Mapping

Each (module × params) combination in a gather stage becomes one Snakemake rule with all provider outputs enumerated as named inputs:

```python
rule summary_S1_default:
    input:
        input_0="data/D1/.f4a3b2c1/methods_fast/M1/.29b6dbbe/D1_M1_result.json",
        input_1="data/D2/.a1b2c3d4/methods_fast/M1/.29b6dbbe/D2_M1_result.json",
        input_2="data/D1/.f4a3b2c1/methods_accurate/M2/.7e8f9a0b/D1_M2_result.json",
        input_3="data/D2/.a1b2c3d4/methods_accurate/M2/.7e8f9a0b/D2_M2_result.json",
    output:
        "summary/S1/.default/report.json"
    shell:
        "... --method {input.input_0} {input.input_1} {input.input_2} {input.input_3} ..."
```

The `benchmark:` directive is omitted for gather nodes (same as existing `metric_collectors`).

### 7.7 Output Path Templates

`provides` labels and structural attributes are available as template variables in the `path:` field of `outputs:`.

#### Available variables

| Variable | Resolves to |
|----------|-------------|
| `{label}` | Module ID of the ancestor (or self) that `provides: {label: ...}` |
| `{params.key}` | Value of parameter `key` for the current node |
| `{module.id}` | Module ID of the current node |
| `{module.stage}` | Stage ID of the current node |
| `{module.parent.id}` | Module ID of the input/parent node |
| `{dataset}` | Backward-compat alias; equivalent to `{dataset}` provides label |

#### Variable availability by stage type

| Variable | First map stage | Downstream map stage | Gather stage |
|----------|-----------------|----------------------|--------------|
| `{label}` (provides) | own module ID if self-provides | ancestor module ID | not available |
| `{params.key}` | yes | yes | yes |
| `{module.*}` | yes | yes | `{module.id}`, `{module.stage}` only |
| `{dataset}` (legacy) | own module ID | propagated from ancestor | not available |

#### Examples

```yaml
# Provider stage: {dataset} = own module ID
- id: data
  provides:
    dataset: data.raw
  outputs:
    - id: data.raw
      path: "{dataset}_data.csv"       # → D1_data.csv, D2_data.csv

# Downstream map: {dataset} from ancestor, {method} = own
- id: methods
  provides:
    method: methods.result
  inputs: [data.raw]
  outputs:
    - id: methods.result
      path: "{dataset}_{method}_result.csv"   # → D1_M1_result.csv, etc.

# Gather stage: no provides-derived variables; {params.*} allowed
- id: summary
  inputs:
    - gather: method
  modules:
    - id: reporter
      parameters:
        - format: [html, pdf]
  outputs:
    - id: summary.report
      path: "report.{params.format}"   # → report.html, report.pdf
```

### 7.8 Ordering and Validation

Stages are still processed in document order. Provider stages must appear before their gather consumers. The compiler validates this and errors with a clear message if violated:

- `"No stage provides 'X'"` — `gather: X` references unknown label
- `"Stage 'Y' gathers 'X' but provider stage 'Z' appears after it"` — ordering violation
- `"Gather stage 'X' cannot mix regular and gather inputs"` — constraint violation

### 7.9 Migration from `metric_collectors`

`metric_collectors` continue to work through the existing resolution path. Gather stages are a superset: composable, chainable, and part of the main stage DAG.

Conceptual equivalence:

```yaml
# Old
metric_collectors:
  - id: MC1
    inputs: [methods.result]
    outputs: [{ id: metrics.summary, path: "metrics.json" }]
    repository: ...

# New
stages:
  - id: metrics
    inputs:
      - gather: method
    outputs:
      - id: metrics.summary
        path: "metrics.json"
    modules:
      - id: MC1
        repository: ...
```

### 7.10 Planned Extensions

| Feature | Description | Status |
|---------|-------------|--------|
| `{dataset}` deprecation | Warnings then removal of hardcoded magic | Planned |
| Real toposort | `graphlib.TopologicalSorter` for arbitrary DAG shapes | Planned |
| `expand_output_path()` deprecation | Remove legacy path injection | Planned |
| Partial collapse | `gather: method, over: [method]` — collapse one dim, keep others | Planned |
| Scoped provides | `provides: [{method: output_id}]` — map label to specific output | Planned |
| Mixed map+gather | Same stage maps one dimension, gathers another | Planned |

---

## 8. Resource Allocation (v0.4+)

> **Status:** Implemented.

### 8.1 Overview

Stages, modules, and metric collectors can declare resource requirements via an optional `resources:` block. These are translated directly into Snakemake `resources:` directives, enabling efficient parallel scheduling.

### 8.2 Schema

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

### 8.3 Placement

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

### 8.4 Resolution Hierarchy

Resources are resolved with **most-specific-wins** precedence:

1. Module-level `resources:` (highest priority)
2. Stage-level `resources:`
3. Global default: `cores=2` (applied when no `resources:` block is present anywhere in the chain)

### 8.5 Snakemake Mapping

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

### 8.6 Controlling parallelism at runtime

Pass `--resources` to the generated Snakemake invocation to cap total consumption:

```bash
# Cap to 100 cores and 128 GiB RAM across all concurrent rules
snakemake --cores 100 --resources mem_mb=131072
```

Snakemake's greedy scheduler will run as many rules in parallel as fit within the declared resource pool.

### 8.7 Complete example

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
