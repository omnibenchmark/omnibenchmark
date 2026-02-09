# 004: Omnibenchmark YAML Specification (v0.3)

[![Status: Under Review](https://img.shields.io/badge/Status-Accepted-green.svg)](https://github.com/omnibenchmark/docs/design)
[![Version: 0.4](https://img.shields.io/badge/Version-0.3-blue.svg)](https://github.com/omnibenchmark/docs/design)

**Authors**: ben,
**Date**: 2025-01-20
**Status**: Under Review
**Version**: 0.1
**Supersedes**: N/A
**Reviewed-by**: TBD
**Related Issues**: TBD

## Changes

| Version | Date | Description | Author |
|---------|------|-------------|--------|
| 0.1 | 2026-01-20 | Initial specification for v0.4.0 | ben |

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

## 6. References

1. [Omnibenchmark Documentation](https://docs.omnibenchmark.org)
2. [Snakemake Documentation](https://snakemake.readthedocs.io)
3. [Semantic Versioning](https://semver.org)
