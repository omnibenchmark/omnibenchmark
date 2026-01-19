"""Unit tests for IOFile path validation warnings."""

import pytest
import tempfile
import warnings
from pathlib import Path

from omnibenchmark.model import IOFile, Benchmark


class TestIOFilePathValidation:
    """Test IOFile path validation for old-style full paths."""

    def test_simple_filename_no_warning(self):
        """Test that simple filenames don't trigger warnings."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")

            # Simple filename should not trigger warning
            io_file = IOFile(id="test.file", path="result.txt")

            assert len(w) == 0
            assert io_file.path == "result.txt"

    def test_old_style_full_path_with_slashes_no_warn_without_context(self):
        """Test that paths with slashes don't warn without YAML context."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")

            # When created directly (not from YAML), no context exists so no warning
            io_file = IOFile(
                id="test.file", path="{input}/{stage}/{module}/{params}/result.txt"
            )

            # No warning without context
            assert len(w) == 0
            assert io_file.path == "{input}/{stage}/{module}/{params}/result.txt"

    def test_relative_path_with_slash_no_warn_without_context(self):
        """Test that paths with slashes don't warn without context."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")

            # Any path with slash - no warning without context
            IOFile(id="test.file", path="subdir/result.txt")

            assert len(w) == 0

    def test_filename_with_extension_no_warning(self):
        """Test that filenames with various extensions don't warn."""
        test_paths = [
            "data.csv",
            "results.json",
            "output.txt",
            "metrics.tsv",
            "plot.png",
        ]

        for path in test_paths:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")

                io_file = IOFile(id="test.file", path=path)

                assert len(w) == 0
                assert io_file.path == path

    def test_empty_path_validation(self):
        """Test that empty paths are still caught by validation."""
        with pytest.raises(ValueError):
            IOFile(id="test.file", path="")

    def test_whitespace_only_path_validation(self):
        """Test that whitespace-only paths are caught."""
        with pytest.raises(ValueError):
            IOFile(id="test.file", path="   ")

    def test_metric_collector_paths_with_slashes_no_warning(self):
        """Test that metric collector outputs with slashes don't warn."""
        # When IOFile is created without stage context (e.g., for metric collectors),
        # it should not warn even if path contains slashes
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")

            # This simulates a metric collector output path
            io_file = IOFile(id="metrics.file", path="{input}/metrics/results.json")

            # Should not warn because there's no stage context
            assert len(w) == 0
            assert io_file.path == "{input}/metrics/results.json"


class TestBenchmarkPathValidationIntegration:
    """Integration tests for path validation in full benchmark YAML."""

    def test_benchmark_with_stage_full_paths_warns(self):
        """Test that stage outputs with full paths trigger warnings."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            f.write("""
id: test_benchmark
description: Test benchmark
version: "1.0.0"
benchmarker: "Test User <test@example.com>"
api_version: "0.3.0"
software_backend: host
software_environments:
  python:
    description: "Python 3.12"
    conda: envs/python.yaml
stages:
  - id: data
    modules:
      - id: test_module
        name: "Test Module"
        software_environment: "python"
        repository:
          url: https://github.com/user/repo
          commit: abc123
    outputs:
      - id: data.result
        path: "{input}/{stage}/{module}/{params}/result.txt"
""")
            yaml_path = f.name

        try:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")

                _ = Benchmark.from_yaml(Path(yaml_path))

                # Should trigger warning for stage output with full path
                future_warnings = [
                    warning
                    for warning in w
                    if issubclass(warning.category, FutureWarning)
                ]
                assert len(future_warnings) >= 1
                assert any(
                    "contains '/'" in str(warning.message)
                    for warning in future_warnings
                )
        finally:
            Path(yaml_path).unlink()

    def test_benchmark_with_metric_collector_full_paths_no_warning(self):
        """Test that metric collector outputs with full paths don't trigger warnings."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            f.write("""
id: test_benchmark
description: Test benchmark
version: "1.0.0"
benchmarker: "Test User <test@example.com>"
api_version: "0.3.0"
software_backend: host
software_environments:
  python:
    description: "Python 3.12"
    conda: envs/python.yaml
stages:
  - id: data
    modules:
      - id: test_module
        name: "Test Module"
        software_environment: "python"
        repository:
          url: https://github.com/user/repo
          commit: abc123
    outputs:
      - id: data.result
        path: "result.txt"
metric_collectors:
  - id: metrics
    description: "Metric collector"
    software_environment: "python"
    repository:
      url: https://github.com/user/metrics
      commit: def456
    inputs: ["data.result"]
    outputs:
      - id: metrics.json
        path: "{input}/metrics/metrics.json"
""")
            yaml_path = f.name

        try:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")

                _ = Benchmark.from_yaml(Path(yaml_path))

                # Should NOT trigger warning for metric collector output
                future_warnings = [
                    warning
                    for warning in w
                    if issubclass(warning.category, FutureWarning)
                ]
                path_warnings = [
                    warning
                    for warning in future_warnings
                    if "contains '/'" in str(warning.message)
                ]
                assert len(path_warnings) == 0
        finally:
            Path(yaml_path).unlink()


class TestFactoryPathValidation:
    """Test that factory functions don't produce path validation warnings."""

    def test_make_iofile_default_no_warning(self):
        """Test that make_iofile() doesn't produce warnings by default."""
        from tests.model.factories import make_iofile

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")

            io_file = make_iofile()

            # Should not trigger any warnings with default values
            assert len(w) == 0
            assert io_file.path == "output.txt"
            assert "/" not in io_file.path

    def test_make_stage_default_no_warning(self):
        """Test that make_stage() doesn't produce warnings by default."""
        from tests.model.factories import make_stage

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")

            stage = make_stage()

            # Should not trigger any warnings with default values
            path_warnings = [
                warning for warning in w if "contains '/'" in str(warning.message)
            ]
            assert len(path_warnings) == 0

            # Check that outputs use simple filenames
            for output in stage.outputs:
                assert (
                    "/" not in output.path
                ), f"Output path '{output.path}' should not contain '/'"

    def test_make_benchmark_default_no_warning(self):
        """Test that make_benchmark() doesn't produce warnings by default."""
        from tests.model.factories import make_benchmark

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")

            benchmark = make_benchmark()

            # Should not trigger any path-related warnings with default values
            path_warnings = [
                warning for warning in w if "contains '/'" in str(warning.message)
            ]
            assert len(path_warnings) == 0

            # Check that all stage outputs use simple filenames
            for stage in benchmark.stages:
                for output in stage.outputs:
                    assert (
                        "/" not in output.path
                    ), f"Output path '{output.path}' should not contain '/'"

    def test_factory_explicit_slash_path_produces_warning_in_stage(self):
        """Test that factory functions with explicit slash paths produce warnings in stages."""
        from tests.model.factories import make_benchmark

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")

            # Explicitly create a benchmark with a stage output that has a slash
            _benchmark = make_benchmark(
                stages=[
                    {
                        "id": "test_stage",
                        "modules": [],
                        "outputs": [{"id": "output1", "path": "results/output.txt"}],
                    }
                ]
            )

            # Should trigger warning for stage output with slash
            path_warnings = [
                warning for warning in w if "contains '/'" in str(warning.message)
            ]
            assert len(path_warnings) >= 1
