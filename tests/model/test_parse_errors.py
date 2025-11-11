"""Tests for benchmark parsing error handling and line tracking."""

import warnings

import pytest

from omnibenchmark.model.benchmark import Benchmark, load_yaml_with_lines
from omnibenchmark.model.validation import BenchmarkParseError


@pytest.mark.short
class TestBenchmarkParseError:
    """Tests for BenchmarkParseError exception class."""

    def test_error_creation_with_all_fields(self, tmp_path):
        """Test creating error with all fields."""
        yaml_file = tmp_path / "test.yaml"
        original_error = ValueError("original error")
        values = ["--method", "test", "--threshold", 0.5]

        error = BenchmarkParseError(
            message="Failed to process parameter",
            yaml_file=yaml_file,
            line_number=42,
            stage_id="clustering",
            module_id="genieclust",
            parameter_index=0,
            values=values,
            original_error=original_error,
        )

        assert error.message == "Failed to process parameter"
        assert error.yaml_file == yaml_file
        assert error.line_number == 42
        assert error.stage_id == "clustering"
        assert error.module_id == "genieclust"
        assert error.parameter_index == 0
        assert error.values == values
        assert error.original_error == original_error

    def test_error_can_be_raised_and_caught(self):
        """Test that error can be raised and caught."""
        with pytest.raises(BenchmarkParseError) as excinfo:
            raise BenchmarkParseError(message="Test error")

        assert excinfo.value.message == "Test error"


@pytest.mark.short
class TestYAMLLineTracking:
    """Tests for YAML line number tracking functionality."""

    def test_basic_yaml_loading_with_lines(self, tmp_path):
        """Test loading basic YAML with line numbers."""
        yaml_content = """id: test
description: test description
version: "1.0"
"""
        yaml_file = tmp_path / "test.yaml"
        yaml_file.write_text(yaml_content)

        data, line_map = load_yaml_with_lines(yaml_file)

        # Data should be loaded correctly
        assert data["id"] == "test"
        assert data["description"] == "test description"
        assert data["version"] == "1.0"

        # Line map should exist
        assert isinstance(line_map, dict)

    def test_line_tracking_preserves_data_integrity(self, tmp_path):
        """Test that line tracking doesn't corrupt data."""
        yaml_content = """id: test_benchmark
description: Test description
version: "1.0"
benchmarker: "Test User"
benchmark_yaml_spec: "0.3.0"
software_backend: "conda"
software_environments:
  - id: env1
    conda: env1.yaml
stages:
  - id: stage1
    modules:
      - id: module1
        software_environment: env1
        repository:
          url: https://github.com/test/repo.git
          commit: abc123
    outputs: []
outputs: []
"""
        yaml_file = tmp_path / "test.yaml"
        yaml_file.write_text(yaml_content)

        # Load with line tracking
        data, line_map = load_yaml_with_lines(yaml_file)

        # Verify data integrity
        assert data["id"] == "test_benchmark"
        assert data["description"] == "Test description"
        assert data["version"] == "1.0"
        assert data["benchmarker"] == "Test User"
        assert len(data["software_environments"]) == 1
        assert len(data["stages"]) == 1
        assert len(data["stages"][0]["modules"]) == 1


@pytest.mark.short
class TestParameterParsingErrors:
    """Tests for error handling in parameter parsing."""

    def test_mixed_type_parameters_work(self, tmp_path):
        """Test that mixed type parameters are handled correctly (the fix for #188)."""
        yaml_content = """id: test_benchmark
description: Test benchmark
version: "1.0"
benchmarker: "Test"
benchmark_yaml_spec: "0.3.0"
software_backend: "conda"
software_environments:
  - id: test_env
    conda: test.yaml
stages:
  - id: clustering
    modules:
      - id: genieclust
        software_environment: test_env
        repository:
          url: https://github.com/test/repo.git
          commit: abc123
        parameters:
          - values: ["--method", "genie", "--threshold", 0.5, "--count", 10]
    outputs: []
outputs: []
"""
        yaml_file = tmp_path / "test.yaml"
        yaml_file.write_text(yaml_content)

        # Should not raise an error - mixed types should be handled
        benchmark = Benchmark.from_yaml(yaml_file)

        # Verify parameter was parsed and ID was generated
        assert len(benchmark.stages[0].modules[0].parameters) == 1
        param = benchmark.stages[0].modules[0].parameters[0]
        assert param.id is not None
        assert len(param.id) == 8  # SHA256 hash truncated to 8 chars

    def test_parameter_parsing_error_includes_context(self, tmp_path):
        """Test that parameter parsing errors include helpful context information."""
        import unittest.mock as mock

        yaml_content = """id: test_benchmark
description: Test
version: "1.0"
benchmarker: "Test"
benchmark_yaml_spec: "0.3.0"
software_backend: "conda"
software_environments:
  - id: test_env
    conda: test.yaml
stages:
  - id: clustering
    modules:
      - id: genieclust
        software_environment: test_env
        repository:
          url: https://github.com/test/repo.git
          commit: abc123
        parameters:
          - values: ["--method", "genie", "--threshold", 0.5]
    outputs: []
outputs: []
"""
        yaml_file = tmp_path / "test.yaml"
        yaml_file.write_text(yaml_content)

        # Simulate an error during parameter processing by mocking sorted()
        def mock_sorted(items):
            if any(isinstance(x, str) and x.startswith("--") for x in items):
                raise TypeError(
                    "'<' not supported between instances of 'float' and 'str'"
                )
            return sorted(items)

        with mock.patch("builtins.sorted", side_effect=mock_sorted):
            with pytest.raises(BenchmarkParseError) as excinfo:
                Benchmark.from_yaml(yaml_file)

            # Verify error has all the context needed
            error = excinfo.value
            assert error.message == "Failed to process parameter values"
            assert error.yaml_file == yaml_file
            assert error.stage_id == "clustering"
            assert error.module_id == "genieclust"
            assert error.parameter_index == 0
            # Values should be converted to strings (conversion happens before error)
            assert error.values == ["--method", "genie", "--threshold", "0.5"]
            assert error.original_error is not None
            # Line number should be tracked (we don't care about actual value)
            assert error.line_number is not None


@pytest.mark.short
class TestTopLevelFieldValidationErrors:
    """Tests for error handling in top-level field validation."""

    def test_version_field_validation_error_includes_line_context(self, tmp_path):
        """Test that validation errors for top-level fields include line numbers."""
        yaml_content = """id: test_benchmark
description: Test benchmark with invalid version format
version: "not-a-semver"
benchmarker: "Test User"
benchmark_yaml_spec: "0.3.0"
software_backend: "conda"
software_environments:
  - id: test_env
    conda: test.yaml
stages: []
outputs: []
"""
        yaml_file = tmp_path / "test.yaml"
        yaml_file.write_text(yaml_content)

        # This should raise a BenchmarkParseError due to invalid version format
        with pytest.raises(BenchmarkParseError) as excinfo:
            Benchmark.from_yaml(yaml_file)

        # Verify error has context
        error = excinfo.value
        assert "does not follow strict semantic versioning format" in error.message
        assert error.yaml_file == yaml_file
        assert error.line_number == 3  # version is on line 3
        assert error.original_error is not None

    def test_nested_field_validation_error_includes_context(self, tmp_path):
        """Test that validation errors for nested fields (e.g., storage.endpoint) include context."""
        yaml_content = """id: test_benchmark
description: Test benchmark
version: "1.0"
benchmarker: "Test User"
benchmark_yaml_spec: "0.3.0"
software_backend: "conda"
software_environments:
  - id: test_env
    conda: test.yaml
stages: []
outputs: []
storage:
  endpoint: 12345
"""
        yaml_file = tmp_path / "test.yaml"
        yaml_file.write_text(yaml_content)

        # This should raise a BenchmarkParseError due to endpoint being an int instead of string
        with pytest.raises(BenchmarkParseError) as excinfo:
            Benchmark.from_yaml(yaml_file)

        # Verify error has context
        error = excinfo.value
        assert "Input should be a valid string" in error.message
        assert error.yaml_file == yaml_file
        # Line number should be tracked
        assert error.line_number is not None
        assert error.original_error is not None


@pytest.mark.short
class TestDeprecationWarningForNumericFields:
    """Tests for deprecation warnings when version or benchmark_yaml_spec are numeric."""

    def test_version_as_float_triggers_deprecation_warning(self, tmp_path):
        """Test that numeric version triggers a future warning but still works."""
        yaml_content = """id: test_benchmark
description: Test benchmark with numeric version
version: 1.4
benchmarker: "Test User"
benchmark_yaml_spec: "0.3.0"
software_backend: "conda"
software_environments:
  - id: test_env
    conda: test.yaml
stages: []
outputs: []
"""
        yaml_file = tmp_path / "test.yaml"
        yaml_file.write_text(yaml_content)

        # Should trigger a DeprecationWarning but still load successfully
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            benchmark = Benchmark.from_yaml(yaml_file)

            # Verify warning was triggered
            assert len(w) == 1
            assert issubclass(w[0].category, FutureWarning)
            assert "version" in str(w[0].message)
            assert "should be a string" in str(w[0].message)
            assert "1.4" in str(w[0].message)
            assert "will not be valid in a future release" in str(w[0].message)

        # Verify the benchmark loaded with the converted value
        assert benchmark.version == "1.4"
        assert benchmark.id == "test_benchmark"

    def test_version_as_int_triggers_deprecation_warning(self, tmp_path):
        """Test that integer version triggers a deprecation warning but fails validation."""
        yaml_content = """id: test_benchmark
description: Test benchmark with integer version
version: 2
benchmarker: "Test User"
benchmark_yaml_spec: "0.3.0"
software_backend: "conda"
software_environments:
  - id: test_env
    conda: test.yaml
stages: []
outputs: []
"""
        yaml_file = tmp_path / "test.yaml"
        yaml_file.write_text(yaml_content)

        # Should trigger a DeprecationWarning and then fail because "2" is not valid semver
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")

            with pytest.raises(BenchmarkParseError) as exc_info:
                Benchmark.from_yaml(yaml_file)

            # Verify deprecation warning was triggered before validation error
            assert len(w) == 1
            assert issubclass(w[0].category, FutureWarning)
            assert "version" in str(w[0].message)

            # Verify the error is about semantic versioning format
            assert "does not follow strict semantic versioning format" in str(
                exc_info.value
            )

    def test_benchmark_yaml_spec_as_float_triggers_deprecation_warning(self, tmp_path):
        """Test that numeric benchmark_yaml_spec triggers a deprecation warning but still works."""
        yaml_content = """id: test_benchmark
description: Test benchmark
version: "1.0"
benchmarker: "Test User"
benchmark_yaml_spec: 0.3
software_backend: "conda"
software_environments:
  - id: test_env
    conda: test.yaml
stages: []
outputs: []
"""
        yaml_file = tmp_path / "test.yaml"
        yaml_file.write_text(yaml_content)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            benchmark = Benchmark.from_yaml(yaml_file)

            assert len(w) == 1
            assert issubclass(w[0].category, FutureWarning)
            assert "benchmark_yaml_spec" in str(w[0].message)
            assert "should be a string" in str(w[0].message)

        assert benchmark.benchmark_yaml_spec == "0.3"

    def test_version_as_string_no_warning(self, tmp_path):
        """Test that string version does not trigger a warning."""
        yaml_content = """id: test_benchmark
description: Test benchmark with proper string version
version: "1.4"
benchmarker: "Test User"
benchmark_yaml_spec: "0.3.0"
software_backend: "conda"
software_environments:
  - id: test_env
    conda: test.yaml
stages: []
outputs: []
"""
        yaml_file = tmp_path / "test.yaml"
        yaml_file.write_text(yaml_content)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            benchmark = Benchmark.from_yaml(yaml_file)

            # No future warnings should be triggered
            future_warnings = [x for x in w if issubclass(x.category, FutureWarning)]
            assert len(future_warnings) == 0

        assert benchmark.version == "1.4"
