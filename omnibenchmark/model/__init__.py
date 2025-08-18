"""Pydantic models for Omnibenchmark."""

from omnibenchmark.model.benchmark import (
    # Base classes
    IdentifiableEntity,
    DescribableEntity,
    # Enums
    APIVersion,
    SoftwareBackendEnum,
    RepositoryType,
    StorageAPIEnum,
    # Core models
    Repository,
    Storage,
    Parameter,
    IOFile,
    InputCollection,
    SoftwareEnvironment,
    SoftwareEnvironmentReference,
    Module,
    MetricCollector,
    Stage,
    Benchmark,
    # Exceptions (moved to validation module)
    # ValidationError,
    # Utility functions
    expand_output_path,
    validate_non_empty_string,
    validate_non_empty_commit,
    validate_hex_string,
)

# Converter classes removed - use Benchmark class directly instead
from omnibenchmark.model.validation import (
    ValidationError,
    BenchmarkValidator,
)
from omnibenchmark.model.module import (
    DerivedSoftware,
    ModuleMetadata,
)

__all__ = [
    # Base classes
    "IdentifiableEntity",
    "DescribableEntity",
    # Enums
    "APIVersion",
    "SoftwareBackendEnum",
    "RepositoryType",
    "StorageAPIEnum",
    # Core models
    "Repository",
    "Storage",
    "Parameter",
    "IOFile",
    "InputCollection",
    "SoftwareEnvironment",
    "SoftwareEnvironmentReference",
    "Module",
    "MetricCollector",
    "Stage",
    "Benchmark",
    # Exceptions
    "ValidationError",
    # Validators
    "BenchmarkValidator",
    # Module metadata
    "DerivedSoftware",
    "ModuleMetadata",
    # Utility functions
    "expand_output_path",
    "validate_non_empty_string",
    "validate_non_empty_commit",
    "validate_hex_string",
]
