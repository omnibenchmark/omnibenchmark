"""
Benchmark definition — the typed data model and its validation.

Represents: what a benchmark *is*, as passive Pydantic structures (stages,
modules, parameters, storage config, ...) plus the rules that validate them.
Layer: foundation (data model)
Depends on: nothing internal — the base layer everything else builds on.
"""

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
    Provenance,
    Parameter,
    IOFile,
    InputCollection,
    Resources,
    SoftwareEnvironment,
    SoftwareEnvironmentReference,
    Module,
    MetricCollector,
    Stage,
    Benchmark,
    # Utility functions
    expand_output_path,
    validate_non_empty_string,
    validate_non_empty_commit,
    validate_hex_string,
)

from omnibenchmark.model.validation import (
    ValidationError,
    BenchmarkValidator,
)
from omnibenchmark.model.module import (
    DerivedSoftware,
    ModuleMetadata,
)
from omnibenchmark.model.params import Params
from omnibenchmark.model.resolved import (
    TemplateContext,
    ResolvedModule,
    ResolvedNode,
    ResolvedMetricCollector,
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
    "Provenance",
    "Parameter",
    "IOFile",
    "InputCollection",
    "Resources",
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
    # Core data structures
    "Params",
    # Resolved entities
    "TemplateContext",
    "ResolvedModule",
    "ResolvedNode",
    "ResolvedMetricCollector",
    # Utility functions
    "expand_output_path",
    "validate_non_empty_string",
    "validate_non_empty_commit",
    "validate_hex_string",
]
