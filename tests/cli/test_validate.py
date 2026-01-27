"""Tests for the validate CLI commands."""

import pytest
import tempfile
import shutil
from pathlib import Path

from tests.cli.cli_setup import OmniCLISetup


@pytest.fixture
def cli_setup():
    """Set up CLI testing environment."""
    return OmniCLISetup()


@pytest.fixture
def temp_dir():
    """Create a temporary directory for testing."""
    temp_dir = Path(tempfile.mkdtemp())
    yield temp_dir
    shutil.rmtree(temp_dir)


@pytest.fixture
def mock_module_dir(temp_dir):
    """Create a mock module directory with all required files."""
    module_dir = temp_dir / "test_module"
    module_dir.mkdir()

    # Create CITATION.cff
    citation_content = """cff-version: 1.2.0
message: If you use this software, please cite it as below.
title: Test Module
authors:
  - family-names: Doe
    given-names: John
license: MIT
"""
    (module_dir / "CITATION.cff").write_text(citation_content)

    # Create LICENSE
    license_content = """MIT License

Copyright (c) 2023 John Doe

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction.
"""
    (module_dir / "LICENSE").write_text(license_content)

    # Create omnibenchmark.yaml
    omnibenchmark_content = """# OmniBenchmark Module Configuration

entrypoints:
  default: run.sh
"""
    (module_dir / "omnibenchmark.yaml").write_text(omnibenchmark_content)

    return module_dir


@pytest.fixture
def mock_benchmark_yaml(temp_dir):
    """Create a mock benchmark YAML file."""
    benchmark_content = """name: test_benchmark
version: 1.0.0
description: A test benchmark

stages:
  - name: data
    modules:
      - repo: https://github.com/test/data.git
"""
    benchmark_file = temp_dir / "benchmark.yaml"
    benchmark_file.write_text(benchmark_content)
    return benchmark_file


class TestValidatePlanCLI:
    """Test suite for the validate plan command."""

    @pytest.mark.short
    def test_validate_plan_help(self, cli_setup):
        """Test that validate plan command shows help information."""
        result = cli_setup.call(["validate", "plan", "--help"])

        assert result.returncode == 0
        assert "Validate benchmark YAML plan structure" in result.stdout
        assert "BENCHMARK" in result.stdout

    @pytest.mark.short
    def test_validate_plan_nonexistent_file(self, cli_setup, temp_dir):
        """Test validate plan command with non-existent file."""
        result = cli_setup.call(
            ["validate", "plan", "/non/existent/benchmark.yaml"],
            cwd=str(temp_dir),
        )

        assert result.returncode != 0
        assert "does not exist" in result.stderr

    @pytest.mark.short
    def test_validate_plan_not_yaml_file(self, cli_setup, temp_dir):
        """Test validate plan command with non-YAML file."""
        # Create a non-YAML file
        non_yaml_file = temp_dir / "test.txt"
        non_yaml_file.write_text("This is not a YAML file")

        result = cli_setup.call(
            ["validate", "plan", str(non_yaml_file)],
            cwd=str(temp_dir),
        )

        assert result.returncode != 0
        assert "Invalid benchmark input" in result.stderr or "YAML" in result.stderr

    @pytest.mark.short
    def test_validate_plan_invalid_yaml(self, cli_setup, temp_dir):
        """Test validate plan command with invalid YAML syntax."""
        invalid_yaml = temp_dir / "invalid.yaml"
        invalid_yaml.write_text("invalid: yaml: syntax: [")

        result = cli_setup.call(
            ["validate", "plan", str(invalid_yaml)],
            cwd=str(temp_dir),
        )

        assert result.returncode != 0

    @pytest.mark.short
    def test_validate_plan_valid_yaml(self, cli_setup, mock_benchmark_yaml):
        """Test validate plan command with valid YAML."""
        result = cli_setup.call(
            ["validate", "plan", str(mock_benchmark_yaml)],
            cwd=str(mock_benchmark_yaml.parent),
        )

        # The command should execute without errors
        # It may fail for other reasons (like missing modules), but YAML parsing should work
        assert "YAML file format error" not in result.stderr


class TestValidateModuleCLI:
    """Test suite for the validate module command."""

    @pytest.mark.short
    def test_validate_module_help(self, cli_setup):
        """Test that validate module command shows help information."""
        result = cli_setup.call(["validate", "module", "--help"])

        assert result.returncode == 0
        assert "Validate module" in result.stdout
        assert "--strict" in result.stdout
        assert "PATH" in result.stdout

    @pytest.mark.short
    def test_validate_module_nonexistent_path(self, cli_setup):
        """Test validate module command with non-existent path."""
        result = cli_setup.call(["validate", "module", "/non/existent/path"])

        assert result.returncode != 0
        assert "does not exist" in result.stderr or "not found" in result.stderr.lower()

    @pytest.mark.short
    def test_validate_module_path_is_file_not_directory(self, cli_setup, temp_dir):
        """Test validate module command when path is a file instead of directory."""
        test_file = temp_dir / "test.txt"
        test_file.write_text("not a directory")

        result = cli_setup.call(
            ["validate", "module", str(test_file)],
            cwd=str(temp_dir),
        )

        assert result.returncode != 0
        assert "not a directory" in result.stderr.lower() or result.returncode != 0

    @pytest.mark.short
    def test_validate_module_missing_omnibenchmark_yaml(self, cli_setup, temp_dir):
        """Test validate module should fail if there's no omnibenchmark.yaml file."""
        module_dir = temp_dir / "test_module"
        module_dir.mkdir()

        # Create CITATION.cff
        citation_content = """cff-version: 1.2.0
message: If you use this software, please cite it as below.
title: Test Module
authors:
  - family-names: Doe
    given-names: John
license: MIT
"""
        (module_dir / "CITATION.cff").write_text(citation_content)

        # Create LICENSE
        (module_dir / "LICENSE").write_text(
            "MIT License\n\nPermission is hereby granted..."
        )

        # Do NOT create omnibenchmark.yaml

        result = cli_setup.call(
            ["validate", "module", str(module_dir), "--strict"],
            cwd=str(temp_dir),
        )

        # Should fail in strict mode
        assert result.returncode != 0

    @pytest.mark.short
    def test_validate_module_missing_omnibenchmark_yaml_non_strict(
        self, cli_setup, temp_dir
    ):
        """Test validate module should fail even in non-strict mode if omnibenchmark.yaml is missing."""
        module_dir = temp_dir / "test_module"
        module_dir.mkdir()

        # Create CITATION.cff
        citation_content = """cff-version: 1.2.0
message: If you use this software, please cite it as below.
title: Test Module
authors:
  - family-names: Doe
    given-names: John
license: MIT
"""
        (module_dir / "CITATION.cff").write_text(citation_content)

        # Create LICENSE
        (module_dir / "LICENSE").write_text(
            "MIT License\n\nPermission is hereby granted..."
        )

        # Do NOT create omnibenchmark.yaml

        result = cli_setup.call(
            ["validate", "module", str(module_dir)],  # No --strict flag
            cwd=str(temp_dir),
        )

        # Missing omnibenchmark.yaml is a hard error even in non-strict mode
        assert result.returncode != 0
        assert "omnibenchmark.yaml" in result.stdout

    @pytest.mark.short
    def test_validate_module_with_legacy_config_cfg(self, cli_setup, temp_dir):
        """Test validate module with legacy config.cfg (should pass with deprecation warning)."""
        module_dir = temp_dir / "test_module"
        module_dir.mkdir()

        # Create CITATION.cff
        citation_content = """cff-version: 1.2.0
message: If you use this software, please cite it as below.
title: Test Module
authors:
  - family-names: Doe
    given-names: John
license: MIT
"""
        (module_dir / "CITATION.cff").write_text(citation_content)

        # Create LICENSE
        (module_dir / "LICENSE").write_text(
            "MIT License\n\nPermission is hereby granted..."
        )

        # Create legacy config.cfg instead of omnibenchmark.yaml
        config_content = """[MODULE]
name = test_module
version = 1.0.0
"""
        (module_dir / "config.cfg").write_text(config_content)

        result = cli_setup.call(
            ["validate", "module", str(module_dir)],
            cwd=str(temp_dir),
        )

        # Should pass (not fail) but with a deprecation warning
        # In non-strict mode, it should succeed
        assert (
            result.returncode == 0
            or "deprecated" in result.stdout.lower()
            or "config.cfg" in result.stdout.lower()
        )

    @pytest.mark.short
    def test_validate_module_all_valid(self, cli_setup, mock_module_dir, temp_dir):
        """Test validate module command with all valid files."""
        result = cli_setup.call(
            ["validate", "module", str(mock_module_dir)],
            cwd=str(temp_dir),
        )

        # Should pass validation
        assert result.returncode == 0
        assert "passed" in result.stdout.lower() or "✅" in result.stdout

    @pytest.mark.short
    def test_validate_module_strict_mode(self, cli_setup, temp_dir):
        """Test validate module with --strict flag."""
        module_dir = temp_dir / "test_module"
        module_dir.mkdir()

        # Create minimal files with issues
        citation_content = """cff-version: 1.2.0
message: Test message
title: Test Module
# Missing authors - this is an error
license: MIT
"""
        (module_dir / "CITATION.cff").write_text(citation_content)

        result = cli_setup.call(
            ["validate", "module", str(module_dir), "--strict"],
            cwd=str(temp_dir),
        )

        # Should fail in strict mode due to missing authors
        assert result.returncode != 0

    @pytest.mark.short
    def test_validate_module_non_strict_mode(self, cli_setup, temp_dir):
        """Test validate module without --strict flag (warnings only)."""
        module_dir = temp_dir / "test_module"
        module_dir.mkdir()

        # Create minimal files with issues
        citation_content = """cff-version: 1.2.0
message: Test message
title: Test Module
# Missing authors - treated as warning in non-strict mode
license: MIT
"""
        (module_dir / "CITATION.cff").write_text(citation_content)

        # Create LICENSE
        (module_dir / "LICENSE").write_text("MIT License")

        # Create omnibenchmark.yaml
        omnibenchmark_content = """# OmniBenchmark Module Configuration

entrypoints:
  default: run.sh
"""
        (module_dir / "omnibenchmark.yaml").write_text(omnibenchmark_content)

        result = cli_setup.call(
            ["validate", "module", str(module_dir)],
            cwd=str(temp_dir),
        )

        # Should pass in non-strict mode (exit code 0) but show warnings for missing authors
        assert result.returncode == 0
        assert "authors" in result.stdout.lower()

    @pytest.mark.short
    def test_validate_module_missing_citation(self, cli_setup, temp_dir):
        """Test validate module with missing CITATION.cff."""
        module_dir = temp_dir / "test_module"
        module_dir.mkdir()

        # Only create LICENSE and omnibenchmark.yaml
        (module_dir / "LICENSE").write_text("MIT License")
        (module_dir / "omnibenchmark.yaml").write_text("name: test\nversion: 1.0.0")

        result = cli_setup.call(
            ["validate", "module", str(module_dir), "--strict"],
            cwd=str(temp_dir),
        )

        # Should fail due to missing CITATION.cff
        assert result.returncode != 0

    @pytest.mark.short
    def test_validate_module_missing_license(self, cli_setup, temp_dir):
        """Test validate module with missing LICENSE file."""
        module_dir = temp_dir / "test_module"
        module_dir.mkdir()

        # Create CITATION.cff and omnibenchmark.yaml but not LICENSE
        citation_content = """cff-version: 1.2.0
message: Test message
title: Test Module
authors:
  - family-names: Doe
    given-names: John
"""
        (module_dir / "CITATION.cff").write_text(citation_content)
        (module_dir / "omnibenchmark.yaml").write_text("name: test\nversion: 1.0.0")

        result = cli_setup.call(
            ["validate", "module", str(module_dir)],
            cwd=str(temp_dir),
        )

        # Should show warning about missing LICENSE file
        assert "LICENSE" in result.stdout or "license" in result.stdout.lower()

    @pytest.mark.short
    def test_validate_module_invalid_citation_yaml(self, cli_setup, temp_dir):
        """Test validate module with invalid CITATION.cff YAML syntax."""
        module_dir = temp_dir / "test_module"
        module_dir.mkdir()

        # Create invalid YAML
        (module_dir / "CITATION.cff").write_text("invalid: yaml: [")
        (module_dir / "LICENSE").write_text("MIT License")
        (module_dir / "omnibenchmark.yaml").write_text("name: test\nversion: 1.0.0")

        result = cli_setup.call(
            ["validate", "module", str(module_dir), "--strict"],
            cwd=str(temp_dir),
        )

        # Should fail due to invalid YAML
        assert result.returncode != 0
