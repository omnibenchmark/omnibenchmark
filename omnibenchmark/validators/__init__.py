"""Runtime artifact validation for omnibenchmark module outputs.

This module validates artifacts (output files) produced by benchmark modules at runtime,
complementing other validation layers in omnibenchmark.

Validation Layers in Omnibenchmark:
------------------------------------
1. **Benchmark YAML validation** (omnibenchmark.model.validation.BenchmarkValidator):
   Validates benchmark specification structure, IDs, references, and schema

2. **Module metadata validation** (omnibenchmark.model.module.ModuleMetadata):
   Validates module metadata YAML (name, author, license, derivatives)

3. **Artifact validation** (this module):
   Validates actual output files produced by module execution match expected schemas

Example validation.yaml (placed in module repository):
-------------------------------------------------------
    results.csv:
      - is_file: true
      - has_column: "score"
      - not_na: "score"

    metadata.json:
      - is_file: true
      - has_key: "timestamp"

    model_output.txt.gz:
      - is_file: true
      - checksum: "sha256:abc123..."

Supported Validators:
---------------------
    - is_file: Check if artifact is a file
    - is_dir: Check if artifact is a directory
    - matches: Regex pattern matching on content
    - checksum: File integrity (format: "algorithm:hash")
    - has_column: Check CSV has column header
    - not_na: Check CSV column has no N/A values
    - has_key: Check JSON has top-level key

See Also:
---------
- tests/data/Benchmark_001.yaml: Example benchmark specification
- tests/data/Benchmark_002.yaml: Example benchmark specification
- omnibenchmark.model.validation: Syntactic benchmark validation
"""

from omnibenchmark.validators.engine import (
    # Core orchestration
    run_engine,
    run_validation_for_module,
    apply_rule,
    load_validation_rules,
    # Result classes
    ValidationResult,
    ModuleValidationResult,
    # Exception hierarchy
    ValidationError,
    ArtifactNotFoundError,
    RuleTypeError,
    FileSystemError,
    ContentError,
    ParsingError,
    # Type definitions
    ValidationRule,
    ValidationFileContent,
    # Registry
    VALIDATOR_REGISTRY,
)

__all__ = [
    # Core functions
    "run_engine",
    "run_validation_for_module",
    "apply_rule",
    "load_validation_rules",
    # Result classes
    "ValidationResult",
    "ModuleValidationResult",
    # Exceptions
    "ValidationError",
    "ArtifactNotFoundError",
    "RuleTypeError",
    "FileSystemError",
    "ContentError",
    "ParsingError",
    # Types
    "ValidationRule",
    "ValidationFileContent",
    # Registry
    "VALIDATOR_REGISTRY",
]

__author__ = "ben"
