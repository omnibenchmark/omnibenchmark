"""Tests for validation consistency using centralized factories.

These tests ensure that the Pydantic-based validation properly catches
invalid configurations at object creation time, following the spec.
"""

import pytest
from pydantic import ValidationError

from omnibenchmark.model.validation import (
    ValidationError as OmnibenchmarkValidationError,
)
from .factories import (
    make_benchmark,
    make_minimal_benchmark,
    make_complete_benchmark,
    make_invalid_benchmark_for_testing,
    make_iofile,
)


@pytest.mark.short
class TestValidationConsistency:
    """Test validation behavior with centralized factories."""

    def test_validate_duplicate_stage_ids(self):
        """Test that duplicate stage IDs are caught during object creation."""
        invalid_data = make_invalid_benchmark_for_testing("duplicate_stage_ids")

        with pytest.raises((ValidationError, OmnibenchmarkValidationError)) as exc_info:
            make_benchmark(**invalid_data)

        error_msg = str(exc_info.value)
        assert "duplicate stage ids" in error_msg.lower()

    def test_validate_undefined_environment(self):
        """Test that undefined environment references are caught."""
        invalid_data = make_invalid_benchmark_for_testing("undefined_environment")

        with pytest.raises((ValidationError, OmnibenchmarkValidationError)) as exc_info:
            make_benchmark(**invalid_data)

        error_msg = str(exc_info.value)
        assert "undefined_env" in error_msg

    def test_validate_invalid_input_reference(self):
        """Test that invalid input references are caught."""
        invalid_data = make_invalid_benchmark_for_testing("invalid_input_reference")

        with pytest.raises((ValidationError, OmnibenchmarkValidationError)) as exc_info:
            make_benchmark(**invalid_data)

        error_msg = str(exc_info.value)
        assert "not valid" in error_msg.lower()

    def test_validate_absolute_paths(self):
        """Test that absolute paths are rejected."""
        invalid_data = make_invalid_benchmark_for_testing("absolute_path")

        with pytest.raises((ValidationError, OmnibenchmarkValidationError)) as exc_info:
            make_benchmark(**invalid_data)

        error_msg = str(exc_info.value)
        assert "absolute" in error_msg.lower() or "relative" in error_msg.lower()

    def test_valid_benchmark_creation(self):
        """Test that valid benchmarks can be created without errors."""
        # Test minimal valid benchmark
        minimal = make_minimal_benchmark()
        assert minimal.id == "test_benchmark"
        assert len(minimal.software_environments) >= 1

        # Test complete valid benchmark
        complete = make_complete_benchmark()
        assert complete.id == "test_benchmark"
        assert len(complete.stages) >= 2
        assert len(complete.software_environments) >= 2
        assert complete.metric_collectors is not None
        assert len(complete.metric_collectors) >= 1

    def test_software_environment_validation(self):
        """Test software environment reference validation."""
        # Valid: environment exists
        valid_benchmark = make_benchmark(
            software_environments=[{"id": "test_env", "conda": "env.yaml"}],
            stages=[
                {
                    "id": "stage1",
                    "modules": [{"id": "mod1", "software_environment": "test_env"}],
                    "outputs": [{"id": "out1", "path": "results.txt"}],
                }
            ],
        )
        errors = []
        valid_benchmark._validate_software_environments(errors)
        assert len(errors) == 0

        # Invalid: environment doesn't exist
        with pytest.raises((ValidationError, OmnibenchmarkValidationError)):
            make_benchmark(
                software_environments=[{"id": "existing_env", "conda": "env.yaml"}],
                stages=[
                    {
                        "id": "stage1",
                        "modules": [
                            {"id": "mod1", "software_environment": "missing_env"}
                        ],
                        "outputs": [],
                    }
                ],
            )

    def test_metric_collector_validation(self):
        """Test metric collector validation."""
        # Valid: metric collector references existing environment and valid inputs
        valid_benchmark = make_benchmark(
            software_environments=[{"id": "metrics_env", "conda": "env.yaml"}],
            stages=[
                {
                    "id": "data_stage",
                    "modules": [],
                    "outputs": [{"id": "data_output", "path": "data/output.csv"}],
                }
            ],
            metric_collectors=[
                {
                    "id": "collector1",
                    "software_environment": "metrics_env",
                    "inputs": [{"id": "data_output", "path": "data/output.csv"}],
                    "outputs": [
                        {"id": "metrics_result", "path": "metrics/report.html"}
                    ],
                }
            ],
        )
        assert valid_benchmark.id == "test_benchmark"

        # Invalid: metric collector references non-existent input
        with pytest.raises((ValidationError, OmnibenchmarkValidationError)):
            make_benchmark(
                software_environments=[{"id": "metrics_env", "conda": "env.yaml"}],
                stages=[
                    {
                        "id": "data_stage",
                        "modules": [],
                        "outputs": [{"id": "real_output", "path": "data/output.csv"}],
                    }
                ],
                metric_collectors=[
                    {
                        "id": "collector1",
                        "software_environment": "metrics_env",
                        "inputs": [{"id": "fake_output", "path": "fake/output.csv"}],
                        "outputs": [],
                    }
                ],
            )

    def test_backend_specific_validation(self):
        """Test software backend-specific validation."""
        # Conda backend with conda environment
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

        # Apptainer backend with apptainer environment
        apptainer_benchmark = make_benchmark(
            software_backend="apptainer",
            software_environments=[
                {"id": "container_env", "apptainer": "container.sif"}
            ],
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

        # Host backend (no specific environment required)
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

    def test_file_path_validation(self):
        """Test file path validation rules."""
        # Valid relative paths
        valid_files = [
            make_iofile(path="relative/path/file.txt"),
            make_iofile(path="simple.txt"),
            make_iofile(path="results/deep/nested/output.json"),
        ]
        for file_obj in valid_files:
            assert not file_obj.path.startswith("/")

        # Invalid paths should be caught during benchmark creation
        with pytest.raises((ValidationError, OmnibenchmarkValidationError)):
            make_benchmark(
                stages=[
                    {
                        "id": "bad_stage",
                        "modules": [],
                        "outputs": [{"id": "bad_file", "path": "/absolute/path.txt"}],
                    }
                ]
            )

        with pytest.raises((ValidationError, OmnibenchmarkValidationError)):
            make_benchmark(
                stages=[
                    {
                        "id": "empty_stage",
                        "modules": [],
                        "outputs": [{"id": "empty_file", "path": ""}],
                    }
                ]
            )

    def test_stage_input_output_consistency(self):
        """Test that stage inputs properly reference existing outputs."""
        # Valid: stage input references output from previous stage
        valid_benchmark = make_benchmark(
            stages=[
                {
                    "id": "producer",
                    "modules": [],
                    "outputs": [{"id": "intermediate_data", "path": "temp/data.csv"}],
                },
                {
                    "id": "consumer",
                    "modules": [],
                    "inputs": [["intermediate_data"]],  # InputCollection format
                    "outputs": [{"id": "final_result", "path": "results/final.json"}],
                },
            ]
        )
        assert len(valid_benchmark.stages) == 2

        # Invalid: stage input references non-existent output
        with pytest.raises((ValidationError, OmnibenchmarkValidationError)):
            make_benchmark(
                stages=[
                    {
                        "id": "producer",
                        "modules": [],
                        "outputs": [{"id": "real_output", "path": "temp/data.csv"}],
                    },
                    {
                        "id": "consumer",
                        "modules": [],
                        "inputs": [["fake_output"]],  # References non-existent output
                        "outputs": [],
                    },
                ]
            )

    def test_comprehensive_valid_benchmark(self):
        """Test creation of a comprehensive, fully valid benchmark."""
        benchmark = make_benchmark(
            id="comprehensive_test",
            name="Comprehensive Test Benchmark",
            description="A complete benchmark for testing all features",
            version="2.0.0",
            benchmarker="Test Engineer",
            software_backend="conda",
            software_environments=[
                {"id": "conda_env", "conda": "env.yaml"},
                {"id": "analysis_env", "conda": "analysis.yaml"},
            ],
            stages=[
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
                        {"id": "analysis_mod", "software_environment": "analysis_env"}
                    ],
                    "outputs": [{"id": "results", "path": "results/analysis.json"}],
                },
            ],
            metric_collectors=[
                {
                    "id": "performance_metrics",
                    "software_environment": "conda_env",
                    "inputs": [{"id": "results", "path": "results/analysis.json"}],
                }
            ],
        )

        # Verify all major components are present and valid
        assert benchmark.id == "comprehensive_test"
        assert benchmark.version == "2.0.0"
        assert len(benchmark.software_environments) >= 2
        assert len(benchmark.stages) >= 2
        assert benchmark.metric_collectors is not None
        assert len(benchmark.metric_collectors) >= 1

        # Validate the benchmark structure
        benchmark.validate_model_structure()
        errors = []
        benchmark._validate_software_environments(errors)
        assert len(errors) == 0

        # Verify relationships are consistent
        env_ids = {env.id for env in benchmark.software_environments}
        for stage in benchmark.stages:
            for module in stage.modules:
                assert module.software_environment in env_ids

        for collector in benchmark.metric_collectors or []:
            assert collector.software_environment in env_ids

    def test_error_message_quality(self):
        """Test that validation error messages are helpful and specific."""
        # Test duplicate ID error message
        try:
            invalid_data = make_invalid_benchmark_for_testing("duplicate_stage_ids")
            make_benchmark(**invalid_data)
            pytest.fail("Expected validation error")
        except (ValidationError, OmnibenchmarkValidationError) as e:
            error_msg = str(e)
            assert "duplicate" in error_msg.lower()
            assert "stage" in error_msg.lower()
            assert "duplicate_id" in error_msg

        # Test undefined environment error message
        try:
            invalid_data = make_invalid_benchmark_for_testing("undefined_environment")
            make_benchmark(**invalid_data)
            pytest.fail("Expected validation error")
        except (ValidationError, OmnibenchmarkValidationError) as e:
            error_msg = str(e)
            assert "undefined_env" in error_msg
            assert "not declared" in error_msg.lower()
