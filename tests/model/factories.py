"""Centralized test factories that respect the omnibenchmark YAML spec.

These factories create valid test objects with minimal boilerplate, following
the structure defined in the auto-generated benchmark template.
"""

from typing import Dict, Any
from omnibenchmark.model import (
    Benchmark,
    SoftwareEnvironment,
    SoftwareBackendEnum,
    StorageAPIEnum,
    Repository,
    Module,
    Stage,
    MetricCollector,
    IOFile,
    Storage,
    InputCollection,
)


# =============================================================================
# Core Entity Factories
# =============================================================================


def make_repository(**kwargs) -> Repository:
    """Create a Repository following the spec."""
    defaults = {
        "url": "https://github.com/example/repo.git",
        "commit": "abc123def456",  # Realistic commit hash
        "type": "git",
    }
    return Repository(**{**defaults, **kwargs})


def make_iofile(**kwargs) -> IOFile:
    """Create an IOFile following the spec."""
    defaults = {
        "id": "output1",
        "path": "results/output.txt",  # Always relative paths
    }
    return IOFile(**{**defaults, **kwargs})


def make_storage(**kwargs) -> Storage:
    """Create Storage config following the spec."""
    defaults = {
        "api": StorageAPIEnum.s3,
        "endpoint": "https://storage.example.com",
        "bucket_name": "test-bucket",
    }
    return Storage(**{**defaults, **kwargs})


# =============================================================================
# Software Environment Factories
# =============================================================================


def make_software_environment(**kwargs) -> SoftwareEnvironment:
    """Create a SoftwareEnvironment following the spec."""
    defaults = {
        "id": "python_env",
        "name": "Python Environment",
        "description": "Test Python environment",
        "conda": "environment.yaml",
    }
    return SoftwareEnvironment(**{**defaults, **kwargs})


def make_apptainer_environment(**kwargs) -> SoftwareEnvironment:
    """Create an Apptainer environment."""
    defaults = {
        "id": "apptainer_env",
        "name": "Apptainer Environment",
        "description": "Test Apptainer environment",
        "apptainer": "container.sif",
    }
    return SoftwareEnvironment(**{**defaults, **kwargs})


def make_docker_environment(**kwargs) -> SoftwareEnvironment:
    """Create a Docker environment."""
    defaults = {
        "id": "docker_env",
        "name": "Docker Environment",
        "description": "Test Docker environment",
        "docker": "python:3.9",
    }
    return SoftwareEnvironment(**{**defaults, **kwargs})


# =============================================================================
# Module and Stage Factories
# =============================================================================


def make_module(**kwargs) -> Module:
    """Create a Module following the spec."""
    defaults = {
        "id": "test_module",
        "name": "Test Module",
        "description": "A test module",
        "software_environment": "python_env",
        "repository": make_repository(),
        "outputs": [make_iofile()],
    }

    # Handle nested objects
    if "repository" in kwargs and isinstance(kwargs["repository"], dict):
        kwargs["repository"] = make_repository(**kwargs["repository"])
    if "outputs" in kwargs and kwargs["outputs"]:
        kwargs["outputs"] = [
            make_iofile(**out) if isinstance(out, dict) else out
            for out in kwargs["outputs"]
        ]

    return Module(**{**defaults, **kwargs})


def make_stage(**kwargs) -> Stage:
    """Create a Stage following the spec."""
    defaults = {
        "id": "preprocessing",
        "name": "Data Preprocessing",
        "description": "Preprocessing stage",
        "modules": [make_module()],
        "outputs": [make_iofile()],
    }

    # Handle nested objects
    if "modules" in kwargs and kwargs["modules"] is not None:
        kwargs["modules"] = [
            make_module(**mod) if isinstance(mod, dict) else mod
            for mod in kwargs["modules"]
        ]
    if "outputs" in kwargs and kwargs["outputs"] is not None:
        kwargs["outputs"] = [
            make_iofile(**out) if isinstance(out, dict) else out
            for out in kwargs["outputs"]
        ]
    if "inputs" in kwargs and kwargs["inputs"] is not None:
        kwargs["inputs"] = [
            InputCollection(entries=inp) if isinstance(inp, list) else inp
            for inp in kwargs["inputs"]
        ]

    return Stage(**{**defaults, **kwargs})


def make_metric_collector(**kwargs) -> MetricCollector:
    """Create a MetricCollector following the spec."""
    defaults = {
        "id": "metrics",
        "name": "Performance Metrics",
        "description": "Collect performance metrics",
        "software_environment": "python_env",
        "repository": make_repository(),
        "inputs": [],
        "outputs": [make_iofile(id="metrics_report", path="metrics/report.html")],
    }

    # Handle nested objects
    if "repository" in kwargs and isinstance(kwargs["repository"], dict):
        kwargs["repository"] = make_repository(**kwargs["repository"])
    if "inputs" in kwargs and kwargs["inputs"] is not None:
        kwargs["inputs"] = [
            make_iofile(**inp) if isinstance(inp, dict) else inp
            for inp in kwargs["inputs"]
        ]
    if "outputs" in kwargs and kwargs["outputs"] is not None:
        kwargs["outputs"] = [
            make_iofile(**out) if isinstance(out, dict) else out
            for out in kwargs["outputs"]
        ]

    return MetricCollector(**{**defaults, **kwargs})


# =============================================================================
# Main Benchmark Factory
# =============================================================================


def make_benchmark(**kwargs) -> Benchmark:
    """Create a Benchmark following the spec with smart defaults.

    This factory creates a valid benchmark that respects the schema and
    automatically sets up consistent relationships between entities.

    Example:
        # Minimal benchmark
        benchmark = make_benchmark()

        # Custom benchmark with validation-friendly setup
        benchmark = make_benchmark(
            id="my_benchmark",
            software_backend="conda",
            stages=[{"id": "stage1", "modules": [{"id": "mod1"}]}]
        )
    """
    # Default software environment that will be referenced
    default_env = make_software_environment()

    defaults = {
        "id": "test_benchmark",
        "name": "Test Benchmark",
        "description": "A test benchmark following the spec",
        "benchmarker": "Test User",
        "version": "1.0.0",
        "software_backend": SoftwareBackendEnum.conda,
        "storage": make_storage(),
        "software_environments": [default_env],
        "stages": [make_stage()],
        "metric_collectors": [],
    }

    # Start with defaults, then apply user overrides
    final_kwargs = {**defaults, **kwargs}

    # Handle software_environments - convert to list if needed
    if "software_environments" in kwargs:
        if isinstance(kwargs["software_environments"], list):
            final_kwargs["software_environments"] = [
                make_software_environment(**env) if isinstance(env, dict) else env
                for env in kwargs["software_environments"]
            ]
        elif isinstance(kwargs["software_environments"], dict):
            # Handle dict format from YAML: {"env_id": {...}}
            envs = []
            for env_id, env_config in kwargs["software_environments"].items():
                if isinstance(env_config, dict):
                    env_config = {**env_config, "id": env_id}
                    envs.append(make_software_environment(**env_config))
                else:
                    envs.append(env_config)
            final_kwargs["software_environments"] = envs

    # Get available environment IDs
    _env_ids = [env.id for env in final_kwargs["software_environments"]]

    # Handle stages
    if "stages" in kwargs and kwargs["stages"] is not None:
        final_kwargs["stages"] = [
            make_stage(**stage) if isinstance(stage, dict) else stage
            for stage in kwargs["stages"]
        ]

    # Handle metric_collectors
    if "metric_collectors" in kwargs and kwargs["metric_collectors"] is not None:
        final_kwargs["metric_collectors"] = [
            make_metric_collector(**mc) if isinstance(mc, dict) else mc
            for mc in kwargs["metric_collectors"]
        ]

    return Benchmark(**final_kwargs)


# =============================================================================
# Convenience Functions for Common Test Patterns
# =============================================================================


def make_minimal_benchmark(**kwargs) -> Benchmark:
    """Create a minimal valid benchmark."""
    minimal_defaults = {
        "stages": [],
        "metric_collectors": [],
        "software_environments": [make_software_environment()],
    }
    return make_benchmark(**{**minimal_defaults, **kwargs})


def make_complete_benchmark(**kwargs) -> Benchmark:
    """Create a benchmark with all optional fields populated."""
    complete_defaults = {
        "software_backend": SoftwareBackendEnum.conda,
        "software_environments": [
            make_software_environment(id="conda_env", conda="env.yaml"),
            make_software_environment(id="container_env", conda="container_env.yaml"),
        ],
        "stages": [
            {
                "id": "preprocessing",
                "modules": [
                    {"id": "preprocess_mod", "software_environment": "conda_env"}
                ],
                "outputs": [{"id": "cleaned_data", "path": "data/cleaned.csv"}],
            },
            {
                "id": "analysis",
                "modules": [
                    {"id": "analysis_mod", "software_environment": "container_env"}
                ],
                "inputs": [["cleaned_data"]],
                "outputs": [{"id": "results", "path": "results/analysis.json"}],
            },
        ],
        "metric_collectors": [
            {
                "id": "performance_metrics",
                "software_environment": "conda_env",
                "inputs": [{"id": "results", "path": "results/analysis.json"}],
            }
        ],
    }
    return make_benchmark(**{**complete_defaults, **kwargs})


def make_invalid_benchmark_for_testing(error_type: str, **kwargs) -> Dict[str, Any]:
    """Create benchmark data that will fail validation for testing purposes.

    Returns a dict (not a Benchmark object) so tests can control when validation occurs.
    """
    base_data = {
        "id": "test_benchmark",
        "description": "Test benchmark",
        "benchmarker": "Test User",
        "version": "1.0.0",
        "software_backend": "conda",
        "storage": {"api": "S3", "endpoint": "https://storage.example.com"},
        "software_environments": [{"id": "python_env", "conda": "env.yaml"}],
        "stages": [],
        "metric_collectors": [],
    }

    if error_type == "duplicate_stage_ids":
        base_data["stages"] = [
            {"id": "duplicate_id", "modules": [], "outputs": []},
            {"id": "duplicate_id", "modules": [], "outputs": []},
        ]
    elif error_type == "undefined_environment":
        base_data["stages"] = [
            {
                "id": "stage1",
                "modules": [{"id": "mod1", "software_environment": "undefined_env"}],
                "outputs": [],
            }
        ]
    elif error_type == "invalid_input_reference":
        base_data["metric_collectors"] = [
            {
                "id": "metrics",
                "software_environment": "python_env",
                "inputs": [{"id": "nonexistent_output", "path": "fake.txt"}],
                "outputs": [],
            }
        ]
    elif error_type == "absolute_path":
        base_data["stages"] = [
            {
                "id": "stage1",
                "modules": [],
                "outputs": [{"id": "bad_output", "path": "/absolute/path.txt"}],
            }
        ]

    return {**base_data, **kwargs}


# =============================================================================
# Backwards Compatibility Functions
# =============================================================================


def make_benchmark_dict(**kwargs) -> Dict[str, Any]:
    """Create a benchmark dictionary (not instantiated) for testing."""
    defaults = {
        "id": "test_benchmark",
        "description": "Test benchmark",
        "benchmarker": "Test User",
        "version": "1.0.0",
        "storage": {
            "api": "S3",
            "endpoint": "https://storage.example.com",
            "bucket_name": "test-bucket",
        },
        "software_backend": "conda",
        "software_environments": [{"id": "python_env", "conda": "env.yaml"}],
        "stages": [],
        "metric_collectors": [],
    }
    return {**defaults, **kwargs}


def make_yaml_content(**kwargs) -> str:
    """Create YAML content for a benchmark."""
    import yaml

    data = make_benchmark_dict(**kwargs)
    return yaml.dump(data, default_flow_style=False)
