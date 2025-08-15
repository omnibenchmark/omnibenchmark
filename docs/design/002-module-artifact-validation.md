# 002: Module Artifact Validation

[![Status: Draft](https://img.shields.io/badge/Status-Draft-yellow.svg)](https://github.com/omnibenchmark/docs/design)
[![Version: 0.1](https://img.shields.io/badge/Version-0.1-blue.svg)](https://github.com/omnibenchmark/docs/design)

**Authors**: ben, ...
**Date**: 2025-07-08
**Status**: Draft
**Version**: 0.1
**Supersedes**: N/A
**Reviewed-by**: TBD
**Related Issues**: TBD

## Changes

| Version | Date | Description | Author |
|---------|------|-------------|--------|
| 0.1 | 2025-07-08 | Initial draft | ben |

## 1. Problem Statement

Currently, omnibenchmark validates module metadata (entrypoints, environments, citations) but lacks mechanisms to validate the artifacts produced by modules. This leads to several issues:

- **Silent failures**: Modules may produce empty or malformed outputs due to runtime errors, but workflow execution continues with `--continue-on-error`, touching empty files that appear as "successful" outputs.
- **Format inconsistencies**: Modules within the same stage may produce outputs in different formats or with incompatible schemas, breaking downstream analysis.
- **Data quality issues**: Outputs may contain entirely NA values, incorrect dimensions, or missing required columns, making benchmark results unreliable.
- **Manual debugging**: Benchmark maintainers must manually inspect outputs to identify problematic modules, which is time-consuming and error-prone.

This design addresses the need for a systematic approach to validate module artifacts against expected contracts, ensuring that benchmark results are based on valid, meaningful data.

### Definitions

- **Artifact**: A file or collection of files produced by the execution of a module. These represent the tangible outputs that modules generate during benchmark execution.
- **Stage**: A logical grouping of modules that perform similar functions in the benchmark workflow (e.g., data preprocessing, method execution, metric calculation).
- **Module Identifier**: Each module has a unique assignment to a stage by means of a unique stage and module identifier, ensuring unambiguous identification of artifact producers.

## 2. Design Goals

- **Contract-based validation**: Define explicit input/output contracts for each benchmark stage, establishing clear expectations for module outputs.
- **Extensible validation framework**: Design a minimal but extensible system that can accommodate future validation needs without architectural changes.
- **Post-hoc validation**: Provide validation capabilities that can be run after workflow execution to identify problematic outputs.
- **Automated disqualification**: Enable machine-parseable validation results that can automatically exclude invalid module outputs from benchmark rounds.
- **Minimal initial implementation**: Start with essential validation capabilities and expand incrementally based on real-world needs.

### Non-Goals

- **Runtime validation**: This design does not cover validation during workflow execution or real-time error handling. In principle, this is possible reusing the same mechanism, but it is not the focus of this design.
- **Performance validation**: Computational resource usage or execution time constraint validation is not in scope.
- **Cross-module consistency**: Validation of relationships between outputs from different modules is not addressed. Each module is validated only according to the artifact contract for a given stage, which is defined in a global `validation.yaml` file.
- **Statistical validation**: Statistical distribution of values or more sophisticated quality checks are not included in the initial scope.

## 3. Proposed Solution

We propose a stage-based validation system that defines output contracts at the benchmark level through a `validation.yaml` file. This system introduces the concept of **artifact contracts** - explicit specifications of what constitutes valid output for each stage. Valid output refers to properties that the artifact must satisfy to be considered valid, regardless of the number or format of the file set that the artifact is made of.

### Core Components

1. **Validation Specification File**: A `validation.yaml` file adjacent to the benchmark's `omnibenchmark.yaml` that defines validation rules per stage.
2. **Mini-DSL**: A domain-specific language for expressing common validation patterns.
3. **Format-aware Loaders**: Extensible data loaders that understand different file formats.
4. **Validation Engine**: A module in `omnibenchmark` that is exposed via `ob check`. The validation engine discovers and executes validation rules. User, in this case, is the benchmarker role.

### Implementation Details

#### Validation File Structure

The `validation.yaml` file defines validation rules organized by stage:

```yaml
# validation.yaml
version: 1.0
stages:
  dataset:
    rules:
      - file_pattern: "*.csv"
        validations:
          - type: not_empty
          - type: has_columns
            columns: ["sample_id", "condition", "value"]
          - type: has_shape
            min_rows: 10
            min_cols: 3
      - file_pattern: "*.json"
        validations:
          - type: not_empty

  processing:
    rules:
      - file_pattern: "{module].csv"
        validations:
          - type: not_empty
          - type: has_columns
            columns: ["cluster_id", "probability"]
          - type: not_na
            columns: ["cluster_id"]

  evaluation:
    rules:
      - file_pattern: "{metric].json"
        validations:
          - type: not_empty
          - type: has_shape
            min_rows: 1
            max_rows: 1
```

#### Mini-DSL Specification

The initial DSL supports these validation primitives:

- **`not_empty`**: Ensures the file is not empty (size > 0 bytes)
- **`has_columns`**: For tabular data, validates presence of required columns
- **`has_shape`**: Validates data dimensions (min/max rows/columns)
- **`not_na`**: Ensures specified columns don't contain all NA values

Future extensions may include:

- **`type`**: Validate column data types (int, float, string)
- **`bounds`**: Validate numeric ranges
- **`format`**: Validate string formats (dates, identifiers)

#### Format Detection and Loaders

The system uses file extensions to determine the appropriate data loader:

- **`.csv`**: Pandas CSV loader
- **`.json`**: JSON loader
- **`.tsv`**: Tab-separated values loader
- **Unknown extensions**: Skip validation with warning

#### Validation Command Interface

```bash
# Auto-discover validation.yaml in current directory
omnibenchmark check --from-results out_dir

# Just artifacts validation
omnibenchmark check --from-results out_dir --no-metadata

# Explicit validation file (assumes output in `out`)
omnibenchmark check --rules validation.yaml

# Validate specific stage
omnibenchmark check --stage method --from-results out_dir

# Validate specific module
omnibenchmark check --module kmeans_sklearn --from-results out_dir

# Output machine-readable results
omnibenchmark check --format json
```

#### Validation Results Format

Validation results are structured for both human readability and machine parsing:

```json
{
  "validation_timestamp": "2025-01-27T10:30:00Z",
  "benchmark": "clustbench",
  "total_modules": 15,
  "passed_modules": 12,
  "failed_modules": 3,
  "results": [
    {
      "stage": "method",
      "module": "kmeans_sklearn",
      "status": "passed",
      "validations": [
        {"rule": "not_empty", "status": "passed"},
        {"rule": "has_columns", "status": "passed"}
      ]
    },
    {
      "stage": "method",
      "module": "hierarchical_ward",
      "status": "failed",
      "validations": [
        {"rule": "not_empty", "status": "failed", "error": "File is empty"},
        {"rule": "has_columns", "status": "skipped"}
      ]
    }
  ]
}
```

## 4. Alternatives Considered

### Alternative 1: Module-level Validation
- **Description**: Place validation rules in each module repository.
- **Pros**: Decentralized, module authors control their validation. Some modules may have specific validation requirements.
- **Cons**: Inconsistent validation across modules, harder to enforce benchmark-wide standards, scattered validation logic
- **Reason for rejection**: Validation should be a benchmark-level concern to ensure consistency

### Alternative 2: Runtime Validation
- **Description**: Validate outputs during workflow execution
- **Pros**: Immediate feedback, could prevent downstream errors
- **Cons**: Complex integration with workflow engines, potential performance impact, harder to implement incremental validation
- **Reason for rejection**: Post-hoc validation is simpler and sufficient for current needs at this stage.

### Alternative 3: Schema-based Validation
- **Description**: Use JSON Schema or similar formal schema languages
- **Pros**: Mature tooling, formal specification
- **Cons**: Overkill for simple validation needs, steep learning curve, limited support for domain-specific patterns
- **Reason for rejection**: Mini-DSL provides better balance of simplicity and expressiveness

## 5. Implementation Plan

### Phase 1: Core Infrastructure
- Implement basic validation engine with file discovery
- Support for `not_empty` validation
- CSV format loader
- Command-line interface with auto-discovery

### Phase 2: DSL Expansion
- Add `has_columns` and `has_shape` validations
- Add `not_na` validation
- JSON format loader
- Enhanced error reporting

### Phase 3: Integration & Polish
- Machine-readable output formats
- Integration with existing omnibenchmark tooling
- Documentation and examples
- Performance optimization
- Make validation mandatory when writing a new benchmark yaml

### Testing Strategy

- **Unit tests**: Each validation rule implementation
- **Integration tests**: End-to-end validation scenarios with sample data
- **Validation metrics**:
  - Accuracy of validation rules (no false positives/negatives)
  - Performance on realistic benchmark outputs
  - Usability testing with benchmark maintainers

## 6. Extensibility & Future Work

The design supports several extension points:

- **New validation rules**: Add to DSL by implementing new validation classes
- **New file formats**: Add format-specific loaders. Design plugins to introduce new file formats.
- **Advanced semantics**: Cross-module validation, statistical checks
- **Runtime integration**: Potential future integration with workflow engines
- **Custom validators**: Plugin system for benchmark-specific validation logic

## 7. References

1. [Module Metadata Design Document](001-module-metadata.md)
2. [Pandas DataFrame Validation Libraries](https://pandera.readthedocs.io/)
3. [JSON Schema Specification](https://json-schema.org/)
4. [Snakemake Continue on Error Documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html#continue-on-errors)
