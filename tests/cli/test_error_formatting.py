"""Tests for CLI error formatting functionality."""

import pytest
from omnibenchmark.cli.error_formatting import pretty_print_parse_error
from omnibenchmark.model.validation import BenchmarkParseError


@pytest.mark.short
class TestPrettyPrintParseError:
    """Tests for pretty_print_parse_error function - focused on behavior not formatting details."""

    def test_minimal_error(self):
        """Test formatting error with only message."""
        error = BenchmarkParseError(message="Test error message")
        result = pretty_print_parse_error(error)
        assert "Test error message" in result

    def test_error_with_context_info(self):
        """Test formatting error with context information."""
        error = BenchmarkParseError(
            message="Failed to process parameter",
            stage_id="clustering",
            module_id="genieclust",
            parameter_index=0,
            values=["--method", "genie", "--threshold", 0.5],
            original_error="'<' not supported between instances of 'float' and 'str'",
        )
        result = pretty_print_parse_error(error)

        # Check that key context is included
        assert "Failed to process parameter" in result
        assert "clustering" in result
        assert "genieclust" in result

    def test_error_with_file_location(self, tmp_path):
        """Test formatting error with file location and line number."""
        yaml_content = """id: test_benchmark
description: Test benchmark
version: "1.0"
benchmarker: "Test"
stages:
  - id: preprocessing
    modules:
      - id: preprocess
        parameters:
          - values: ["--method", "test"]
"""
        yaml_file = tmp_path / "test.yaml"
        yaml_file.write_text(yaml_content)

        error = BenchmarkParseError(
            message="Failed to parse parameter",
            yaml_file=yaml_file,
            line_number=10,
        )
        result = pretty_print_parse_error(error)

        # Should include message and file reference
        assert "Failed to parse parameter" in result
        assert str(yaml_file) in result

    def test_error_with_nonexistent_file(self, tmp_path):
        """Test that errors handle nonexistent files gracefully."""
        nonexistent = tmp_path / "nonexistent.yaml"

        error = BenchmarkParseError(
            message="Error message",
            yaml_file=nonexistent,
            line_number=10,
        )
        result = pretty_print_parse_error(error)

        # Should not crash and should include message
        assert "Error message" in result

    def test_error_with_none_values(self):
        """Test that None values in context are handled gracefully."""
        error = BenchmarkParseError(
            message="Error message",
            stage_id=None,
            module_id=None,
            parameter_index=None,
            values=None,
            original_error=None,
        )
        result = pretty_print_parse_error(error)

        # Should not crash
        assert "Error message" in result

    def test_error_with_zero_parameter_index(self):
        """Test that parameter_index=0 is displayed (not treated as falsy)."""
        error = BenchmarkParseError(
            message="Error message",
            parameter_index=0,
        )
        result = pretty_print_parse_error(error)

        assert "Error message" in result
        assert "0" in result  # Parameter index should appear
