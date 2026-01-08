"""Tests for validation error formatting in CLI commands."""

from pathlib import Path

import pytest
from click.testing import CliRunner
from pydantic import ValidationError as PydanticValidationError

from omnibenchmark.cli.run import format_pydantic_errors, run
from omnibenchmark.cli.validate import validate_plan

data = Path(__file__).parent.parent / "data"


@pytest.mark.short
def test_format_pydantic_errors_missing_fields():
    """Test that format_pydantic_errors properly formats missing field errors."""
    # Create a mock Pydantic validation error
    try:
        from omnibenchmark.model import Benchmark

        # This should fail with missing required fields
        Benchmark(**{"id": "test"})
    except PydanticValidationError as e:
        formatted = format_pydantic_errors(e)

        # Check that the formatted error contains helpful information
        assert "Validation failed:" in formatted
        assert "Missing required field:" in formatted
        # Should mention some of the missing fields
        assert any(
            field in formatted
            for field in ["benchmarker", "version", "software_backend", "stages"]
        )


@pytest.mark.short
def test_format_pydantic_errors_invalid_field():
    """Test that format_pydantic_errors properly formats invalid field errors."""
    try:
        from omnibenchmark.model import Benchmark

        # This should fail with invalid version format
        Benchmark(
            **{
                "id": "test",
                "benchmarker": "test",
                "version": "invalid-version",  # Invalid semantic version
                "software_backend": "conda",
                "software_environments": [],
                "stages": [],
            }
        )
    except PydanticValidationError as e:
        formatted = format_pydantic_errors(e)

        # Check that the formatted error contains helpful information
        assert "Validation failed:" in formatted
        assert "version" in formatted.lower()


@pytest.mark.short
def test_validate_plan_with_missing_fields(tmp_path, caplog):
    """Test that validate plan command shows helpful errors for missing required fields."""
    # Create a minimal invalid YAML file
    invalid_yaml = tmp_path / "invalid_benchmark.yaml"
    invalid_yaml.write_text("""
id: test_benchmark
name: Test Benchmark
# Missing required fields: benchmarker, version, software_backend, stages, etc.
""")

    runner = CliRunner()
    result = runner.invoke(validate_plan, [str(invalid_yaml)])

    # Should exit with error
    assert result.exit_code == 1

    # Should show error message in logs (gets converted to BenchmarkParseError)
    log_output = caplog.text
    assert "Error:" in log_output
    assert "Field required" in log_output or "Missing required field" in log_output


@pytest.mark.short
def test_validate_plan_with_invalid_version(tmp_path, caplog):
    """Test that validate plan command shows helpful errors for invalid field values."""
    # Create a YAML file with invalid version
    invalid_yaml = tmp_path / "invalid_version.yaml"
    invalid_yaml.write_text("""
id: test_benchmark
name: Test Benchmark
benchmarker: test@example.com
version: not-a-version
software_backend: conda
software_environments: []
stages: []
""")

    runner = CliRunner()
    result = runner.invoke(validate_plan, [str(invalid_yaml)])

    # Should exit with error
    assert result.exit_code == 1

    # Should show helpful error message about version in logs
    log_output = caplog.text
    assert "Error:" in log_output
    assert "version" in log_output.lower()
    assert "semantic versioning" in log_output.lower()


@pytest.mark.short
def test_run_with_missing_fields(tmp_path, caplog):
    """Test that run command shows helpful errors for missing required fields."""
    # Create a minimal invalid YAML file
    invalid_yaml = tmp_path / "invalid_benchmark.yaml"
    invalid_yaml.write_text("""
id: test_benchmark
name: Test Benchmark
# Missing required fields
""")

    runner = CliRunner()
    result = runner.invoke(run, [str(invalid_yaml)])

    # Should exit with error
    assert result.exit_code == 1

    # Should show helpful error message in logs
    log_output = caplog.text
    assert "Failed to load benchmark:" in log_output
    assert "Field required" in log_output or "Missing required field" in log_output


@pytest.mark.short
def test_format_pydantic_errors_multiple_errors():
    """Test that format_pydantic_errors handles multiple errors correctly."""
    try:
        from omnibenchmark.model import Benchmark

        # This should fail with multiple missing fields
        Benchmark(**{})
    except PydanticValidationError as e:
        formatted = format_pydantic_errors(e)

        # Check that multiple errors are listed
        assert "Validation failed:" in formatted
        assert formatted.count("Missing required field:") > 1

        # Each error should be on a separate line with proper indentation
        lines = formatted.split("\n")
        assert len(lines) > 2  # Header + multiple error lines
        assert all(line.startswith("  - ") for line in lines[1:] if line.strip())
