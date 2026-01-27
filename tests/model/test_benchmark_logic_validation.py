"""Test cases for benchmark logic validation including parsing edge cases."""

import pytest
import warnings
from pydantic import ValidationError
from omnibenchmark.model.validation import (
    ValidationError as OmnibenchmarkValidationError,
)
from tests.model.factories import make_benchmark


class TestBenchmarkLogicValidation:
    """Test benchmark validation logic including parsing edge cases."""

    @pytest.mark.short
    def test_benchmark_yaml_spec_numeric_compatibility(self):
        """Test that benchmark_yaml_spec accepts numeric values with deprecation warnings."""
        # Test float value
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            benchmark = make_benchmark(benchmark_yaml_spec=0.3)

            assert benchmark.benchmark_yaml_spec == "0.3"
            assert len(w) == 1
            assert issubclass(w[0].category, FutureWarning)
            assert "benchmark_yaml_spec" in str(w[0].message)
            assert "should be a string" in str(w[0].message)

        # Test integer value
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            benchmark = make_benchmark(benchmark_yaml_spec=1)

            assert benchmark.benchmark_yaml_spec == "1"
            assert len(w) == 1
            assert issubclass(w[0].category, FutureWarning)
            assert "benchmark_yaml_spec" in str(w[0].message)
            assert "should be a string" in str(w[0].message)

        # Test that strings work without warnings
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            benchmark = make_benchmark(benchmark_yaml_spec="0.3")

            assert benchmark.benchmark_yaml_spec == "0.3"
            assert len(w) == 0  # No warnings for string values

        # Test that None works without warnings
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            benchmark = make_benchmark(benchmark_yaml_spec=None)

            assert benchmark.benchmark_yaml_spec is None
            assert len(w) == 0  # No warnings for None

        # Test that invalid types still raise errors
        with pytest.raises(ValidationError):
            make_benchmark(benchmark_yaml_spec={"version": "0.3"})

    @pytest.mark.short
    def test_metric_collector_string_inputs(self):
        """Test that metric collectors accept string inputs alongside IOFile objects."""
        # Test with string inputs only
        benchmark = make_benchmark(
            stages=[
                {
                    "id": "methods",
                    "modules": [],
                    "outputs": [{"id": "methods.result", "path": "output.json"}],
                }
            ],
            metric_collectors=[
                {
                    "id": "test_collector",
                    "software_environment": "python_env",  # Use default env
                    "repository": {"url": "test.git", "commit": "abc123"},
                    "inputs": ["methods.result"],  # String input
                    "outputs": [],
                }
            ],
        )

        collector = benchmark.metric_collectors[0]
        assert len(collector.inputs) == 1
        assert collector.inputs[0] == "methods.result"
        assert isinstance(collector.inputs[0], str)

        # Test with mixed string and IOFile inputs
        benchmark = make_benchmark(
            stages=[
                {
                    "id": "methods",
                    "modules": [],
                    "outputs": [
                        {"id": "methods.result", "path": "output.json"},
                        {"id": "data.input", "path": "input.json"},
                    ],
                }
            ],
            metric_collectors=[
                {
                    "id": "mixed_collector",
                    "software_environment": "python_env",  # Use default env
                    "repository": {"url": "test.git", "commit": "abc123"},
                    "inputs": [
                        "methods.result",  # String input
                        {"id": "data.input", "path": "input.json"},  # IOFile object
                    ],
                    "outputs": [],
                }
            ],
        )

        collector = benchmark.metric_collectors[0]
        assert len(collector.inputs) == 2
        assert collector.inputs[0] == "methods.result"
        assert isinstance(collector.inputs[0], str)
        assert hasattr(collector.inputs[1], "id")
        assert collector.inputs[1].id == "data.input"
        assert collector.inputs[1].path == "input.json"

        # Test with IOFile inputs only (existing behavior)
        benchmark = make_benchmark(
            stages=[
                {
                    "id": "methods",
                    "modules": [],
                    "outputs": [{"id": "methods.result", "path": "output.json"}],
                }
            ],
            metric_collectors=[
                {
                    "id": "iofile_collector",
                    "software_environment": "python_env",  # Use default env
                    "repository": {"url": "test.git", "commit": "abc123"},
                    "inputs": [{"id": "methods.result", "path": "output.json"}],
                    "outputs": [],
                }
            ],
        )

        collector = benchmark.metric_collectors[0]
        assert len(collector.inputs) == 1
        assert hasattr(collector.inputs[0], "id")
        assert collector.inputs[0].id == "methods.result"
        assert collector.inputs[0].path == "output.json"

    @pytest.mark.short
    def test_metric_collector_validation_with_string_inputs(self):
        """Test that validation works correctly with string inputs in metric collectors."""
        # Test valid string input reference
        benchmark = make_benchmark(
            stages=[
                {
                    "id": "methods",
                    "modules": [],
                    "outputs": [{"id": "methods.result", "path": "output.json"}],
                }
            ],
            metric_collectors=[
                {
                    "id": "valid_collector",
                    "software_environment": "python_env",  # Use default env
                    "repository": {"url": "test.git", "commit": "abc123"},
                    "inputs": ["methods.result"],  # Valid reference
                    "outputs": [],
                }
            ],
        )

        # Should not raise validation errors for valid references
        try:
            benchmark.validate_model_structure()
        except OmnibenchmarkValidationError:
            pytest.fail(
                "Valid string input reference should not raise validation error"
            )

        # Test that the string input is properly stored
        collector = benchmark.metric_collectors[0]
        assert len(collector.inputs) == 1
        assert collector.inputs[0] == "methods.result"
        assert isinstance(collector.inputs[0], str)
