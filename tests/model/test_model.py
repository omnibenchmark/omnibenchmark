"""Tests for the new Pydantic-based model to replace omni-schema."""

import pytest
import yaml

from pydantic import ValidationError

from omnibenchmark.model import (
    # Helper functions
    validate_non_empty_string,
    validate_non_empty_commit,
    validate_hex_string,
    # Enums
    APIVersion,
    SoftwareBackendEnum,
    RepositoryType,
    StorageAPIEnum,
    # Core entities
    IOFile,
    Parameter,
    Benchmark,
)

from .factories import (
    make_benchmark,
    make_benchmark_dict,
    make_yaml_content,
    make_software_environment,
    make_module,
    make_stage,
    make_metric_collector,
    make_iofile,
    make_repository,
)


# Test helper functions
@pytest.mark.short
class TestHelperFunctions:
    def test_validate_non_empty_string(self):
        """Test non-empty string validation."""
        assert validate_non_empty_string("valid") == "valid"

        with pytest.raises(ValueError, match="must be a non-empty string"):
            validate_non_empty_string("")

        with pytest.raises(ValueError, match="must be a non-empty string"):
            validate_non_empty_string("   ")

    def test_validate_non_empty_commit(self):
        """Test commit hash validation."""
        assert validate_non_empty_commit("abc123") == "abc123"

        with pytest.raises(ValueError, match="Commit must be a non-empty string"):
            validate_non_empty_commit("")

    def test_validate_hex_string(self):
        """Test hexadecimal string validation."""
        # Valid hex strings
        assert validate_hex_string("abc123") == "abc123"
        assert validate_hex_string("ABCDEF") == "ABCDEF"
        assert validate_hex_string("0123456789abcdef") == "0123456789abcdef"

        # Invalid hex strings
        with pytest.raises(ValueError, match="must be a valid hexadecimal string"):
            validate_hex_string("xyz")

        with pytest.raises(ValueError, match="must be a valid hexadecimal string"):
            validate_hex_string("abc123g")

        with pytest.raises(ValueError, match="must be a valid hexadecimal string"):
            validate_hex_string("")


# Test enums
@pytest.mark.short
class TestEnums:
    def test_api_version(self):
        """Test APIVersion enum."""
        assert APIVersion.V0_1_0.value == "0.1.0"
        assert APIVersion.V0_2_0.value == "0.2.0"
        assert APIVersion.V0_3_0.value == "0.3.0"

        assert APIVersion.latest() == "0.3.0"
        assert set(APIVersion.supported_versions()) == {"0.1.0", "0.2.0", "0.3.0"}

    def test_software_backend_enum(self):
        """Test SoftwareBackendEnum."""
        assert SoftwareBackendEnum.host.value == "host"
        assert SoftwareBackendEnum.docker.value == "docker"
        assert SoftwareBackendEnum.apptainer.value == "apptainer"
        assert SoftwareBackendEnum.conda.value == "conda"
        assert SoftwareBackendEnum.envmodules.value == "envmodules"

    def test_repository_type(self):
        """Test RepositoryType enum."""
        assert RepositoryType.git.value == "git"

    def test_storage_api_enum(self):
        """Test StorageAPIEnum."""
        assert StorageAPIEnum.s3.value == "S3"


# Test core entities
@pytest.mark.short
class TestCoreEntities:
    def test_repository(self):
        """Test Repository model."""
        repo = make_repository(url="https://github.com/test/repo.git")
        assert repo.url == "https://github.com/test/repo.git"
        assert repo.commit == "abc123def456"  # default from factory

        # Test validation
        with pytest.raises(ValidationError):
            make_repository(url="", commit="abc123")

        with pytest.raises(ValidationError):
            make_repository(url="https://example.com", commit="")

    def test_parameter(self):
        """Test Parameter model."""
        param = Parameter(id="param1", values=["a", "b", "c"])
        assert param.id == "param1"
        assert param.values == ["a", "b", "c"]

        # Test with empty values
        with pytest.raises(ValidationError):
            Parameter(id="param1", values=[])

    def test_io_file(self):
        """Test IOFile model."""
        file = make_iofile(id="output1", path="results/output.csv")
        assert file.id == "output1"
        assert file.path == "results/output.csv"

        # Test validation
        with pytest.raises(ValidationError):
            IOFile(id="output1", path="")

    def test_software_environment_complete(self):
        """Test SoftwareEnvironment with all fields."""
        env = make_software_environment(
            id="complete_env",
            description="Complete environment",
            easyconfig="test.eb",
            envmodule="test/1.0",
            conda="env.yaml",
            apptainer="test.sif",
        )
        assert env.easyconfig == "test.eb"
        assert env.envmodule == "test/1.0"
        assert env.conda == "env.yaml"
        assert env.apptainer == "test.sif"

    def test_module_with_outputs(self):
        """Test Module with outputs."""
        module = make_module(
            id="module_outputs",
            outputs=[
                {"id": "out1", "path": "output1.txt"},
                {"id": "out2", "path": "output2.txt"},
            ],
        )
        assert module.outputs is not None
        assert len(module.outputs) == 2
        assert module.outputs[0].id == "out1"

    def test_has_environment_reference(self):
        """Test has_environment_reference method."""
        module = make_module(software_environment="env1")
        assert module.has_environment_reference("env1") is True
        assert module.has_environment_reference("env2") is False

        collector = make_metric_collector(software_environment="env1")
        assert collector.has_environment_reference("env1") is True
        assert collector.has_environment_reference("env2") is False


# Test Benchmark
@pytest.mark.short
class TestBenchmark:
    def test_minimal_benchmark(self):
        """Test creating a minimal Benchmark."""
        benchmark = make_benchmark()
        assert benchmark.id == "test_benchmark"
        assert benchmark.description == "A test benchmark following the spec"
        assert benchmark.version == "1.0.0"
        assert benchmark.software_backend == SoftwareBackendEnum.conda

    def test_complete_benchmark(self):
        """Test creating a complete Benchmark."""
        benchmark = make_benchmark(
            id="complete_benchmark",
            software_environments={
                "python_env": {
                    "id": "python_env",
                    "description": "Python environment",
                    "conda": "python.yaml",
                },
                "r_env": {
                    "id": "r_env",
                    "description": "R environment",
                    "conda": "r.yaml",
                },
            },
            stages=[
                {
                    "id": "stage1",
                    "modules": [
                        {
                            "id": "mod1",
                            "name": "Module 1",
                            "software_environment": "python_env",
                        }
                    ],
                }
            ],
            metric_collectors=[{"id": "metrics", "software_environment": "python_env"}],
        )
        assert benchmark.id == "complete_benchmark"
        assert len(benchmark.software_environments) == 2
        assert benchmark.metric_collectors is not None
        assert len(benchmark.metric_collectors) == 1
        assert len(benchmark.stages) == 1

    def test_from_yaml_string(self):
        """Test loading Benchmark from YAML string."""
        yaml_content = make_yaml_content(id="yaml_benchmark")
        benchmark = Benchmark.from_yaml(yaml_content)
        assert benchmark.id == "yaml_benchmark"

    def test_from_yaml_file(self, tmp_path):
        """Test loading Benchmark from YAML file."""
        yaml_file = tmp_path / "test.yaml"
        yaml_file.write_text(make_yaml_content(id="file_benchmark"))

        benchmark = Benchmark.from_yaml(yaml_file)
        assert benchmark.id == "file_benchmark"

    def test_from_dict(self):
        """Test creating Benchmark from dictionary."""
        benchmark_dict = make_benchmark_dict(id="dict_benchmark")
        benchmark = Benchmark.from_dict(benchmark_dict)
        assert benchmark.id == "dict_benchmark"

    def test_merge_benchmarks(self):
        """Test merging two benchmarks."""
        benchmark1 = make_benchmark(stages=[{"id": "stage1"}])

        benchmark2 = make_benchmark(
            description="Updated description",
            stages=[{"id": "stage2", "modules": [], "outputs": []}],
            software_environments=[{"id": "test_env", "conda": "test.yaml"}],
            metric_collectors=[{"id": "metrics", "software_environment": "test_env"}],
        )

        merged = benchmark1.merge_with(benchmark2)
        assert merged.description == "Updated description"
        assert len(merged.stages) == 2
        assert merged.metric_collectors is not None
        assert len(merged.metric_collectors) == 1

    def test_upgrade_to_latest(self):
        """Test upgrading benchmark version."""
        old_benchmark = make_benchmark(benchmark_yaml_spec="0.1.0")
        upgraded = old_benchmark.upgrade_to_latest()
        assert upgraded.benchmark_yaml_spec == APIVersion.latest()

    def test_validate_software_environments(self):
        """Test software environment validation."""
        benchmark = make_benchmark(
            software_backend="conda",
            software_environments={
                "env1": {"conda": "env1.yaml"},
                "env2": {"conda": "env2.yaml"},
            },
            stages=[
                {
                    "modules": [
                        {"software_environment": "env1"},
                        {"software_environment": "env2"},
                    ]
                }
            ],
        )

        # Should not raise any errors
        benchmark.validate_software_environments()

        # Test with missing environment reference
        benchmark.stages[0].modules[0].software_environment = "nonexistent_env"

        with pytest.raises(ValueError, match="not defined"):
            benchmark.validate_software_environments()

    def test_validate_unused_environments(self):
        """Test detection of unused software environments."""
        benchmark = make_benchmark(
            software_environments={"used_env": {}, "unused_env": {}},
            stages=[{"modules": [{"software_environment": "used_env"}]}],
        )

        # Should detect unused environment
        with pytest.warns(UserWarning, match="unused_env"):
            benchmark.validate_software_environments()

    def test_software_backend_flexibility(self):
        """Test that different software backends work with appropriate environments."""
        # Test conda backend
        conda_benchmark = make_benchmark(
            software_backend="conda",
            software_environments=[{"id": "conda_env", "conda": "environment.yaml"}],
            stages=[
                {
                    "id": "stage1",
                    "modules": [{"id": "mod1", "software_environment": "conda_env"}],
                    "outputs": [],
                }
            ],
        )
        assert conda_benchmark.software_backend.value == "conda"

        # Test apptainer backend
        apptainer_benchmark = make_benchmark(
            software_backend="apptainer",
            software_environments=[{"id": "container_env", "apptainer": "image.sif"}],
            stages=[
                {
                    "id": "stage1",
                    "modules": [
                        {"id": "mod1", "software_environment": "container_env"}
                    ],
                    "outputs": [],
                }
            ],
        )
        assert apptainer_benchmark.software_backend.value == "apptainer"

        # Test host backend (most flexible)
        host_benchmark = make_benchmark(
            software_backend="host",
            software_environments=[
                {"id": "host_env", "description": "Host environment"}
            ],
            stages=[
                {
                    "id": "stage1",
                    "modules": [{"id": "mod1", "software_environment": "host_env"}],
                    "outputs": [],
                }
            ],
        )
        assert host_benchmark.software_backend.value == "host"


# Test legacy API
@pytest.mark.short
@pytest.mark.skip(reason="Legacy API not implemented")
class TestLegacyAPI:
    def test_all_legacy_methods(self):
        """Test all legacy API methods with deprecation warnings."""
        pass


# Test error cases
@pytest.mark.short
class TestErrorCases:
    def test_invalid_yaml_format(self, tmp_path):
        """Test handling of invalid YAML format."""
        invalid_yaml = tmp_path / "invalid.yaml"
        invalid_yaml.write_text("{ invalid yaml content")

        with pytest.raises(yaml.YAMLError):
            Benchmark.from_yaml(invalid_yaml)

    def test_missing_required_fields(self):
        """Test validation of missing required fields."""
        with pytest.raises(ValidationError) as excinfo:
            Benchmark(
                id="incomplete",
                description="Missing fields",
                name="Incomplete Benchmark",
                software_backend="conda",
                software_environments=[],
                stages=[],
            )

        errors = excinfo.value.errors()
        required_fields = {error["loc"][0] for error in errors}
        assert "version" in required_fields
        assert "benchmarker" in required_fields
        # storage is optional, so don't check for it

    def test_invalid_enum_values(self):
        """Test validation of invalid enum values."""
        with pytest.raises(ValidationError):
            make_benchmark(software_backend="invalid_backend")

    def test_empty_collections(self):
        """Test validation of empty collections."""
        # Empty parameter values
        with pytest.raises(ValidationError):
            Parameter(id="param", values=[])

        # Empty stage modules - should be allowed
        stage = make_stage(modules=[])
        assert stage.modules == []

    def test_empty_file_paths(self):
        """Test that empty paths are rejected."""
        with pytest.raises(ValidationError):
            IOFile(id="file", path="")

        with pytest.raises(ValidationError):
            IOFile(id="file", path="   ")


# Integration tests
@pytest.mark.short
class TestIntegration:
    def test_full_benchmark_workflow(self, tmp_path):
        """Test complete workflow from YAML to validation."""
        # Create YAML content
        yaml_content = """
id: integration_test
description: Integration test benchmark
version: "1.0"
benchmarker: "Integration Tester"
storage:
  api: "S3"
  endpoint: "https://storage.example.com"
  bucket_name: "integration-bucket"
benchmark_yaml_spec: "0.3.0"
software_backend: "conda"
software_environments:
  python_env:
    id: python_env
    description: "Python environment"
    conda: "envs/python.yaml"
stages:
  - id: preprocessing
    modules:
      - id: preprocess_module
        name: "Data Preprocessing"
        software_environment: "python_env"
        repository:
          url: "https://github.com/example/preprocess.git"
          commit: "abc123"
        outputs:
          - id: "cleaned_data"
            path: "output/cleaned.csv"
    outputs:
      - id: "cleaned_data"
        path: "output/cleaned.csv"
metric_collectors:
  - id: performance_metrics
    name: "Performance Metrics"
    software_environment: "python_env"
    repository:
      url: "https://github.com/example/metrics.git"
      commit: "ghi789"
    inputs:
      - id: "cleaned_data"
        path: "output/cleaned.csv"
    outputs:
      - id: "metrics_report"
        path: "metrics/report.html"
outputs:
  - id: "final_report"
    path: "final/report.pdf"
"""
        yaml_file = tmp_path / "integration.yaml"
        yaml_file.write_text(yaml_content)

        # Load and validate
        benchmark = Benchmark.from_yaml(yaml_file)
        assert benchmark.id == "integration_test"

        # Validate software environments
        benchmark.validate_software_environments()
