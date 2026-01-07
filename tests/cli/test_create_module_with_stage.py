"""Tests for create module with --benchmark and --for-stage options."""

import pytest
from click.testing import CliRunner

from omnibenchmark.cli.create import create_module


@pytest.fixture
def sample_benchmark_yaml(tmp_path):
    """Create a sample benchmark YAML file for testing."""
    benchmark_content = """id: TestBenchmark
description: Test benchmark for module creation
version: "1.0"
benchmarker: "Test Suite"
benchmark_yaml_spec: 0.3
software_backend: host

software_environments:
  host:
    description: "Host environment"

stages:
  - id: data
    modules:
      - id: D1
        name: "Dataset 1"
        software_environment: "host"
        repository:
          url: https://github.com/example/data.git
          commit: abc123
    outputs:
      - id: data.raw
        path: "{dataset}_data.json"
      - id: data.meta
        path: "{dataset}_meta.json"

  - id: methods
    modules:
      - id: M1
        name: "Method 1"
        software_environment: "host"
        repository:
          url: https://github.com/example/method.git
          commit: def456
    inputs:
      - entries:
          - data.raw
          - data.meta
    outputs:
      - id: methods.result
        path: "{dataset}_result.json"

  - id: metrics
    modules:
      - id: Met1
        name: "Metric 1"
        software_environment: "host"
        repository:
          url: https://github.com/example/metric.git
          commit: ghi789
    inputs:
      - entries:
          - methods.result
          - data.raw
    outputs:
      - id: metrics.score
        path: "{dataset}_score.json"
"""
    benchmark_file = tmp_path / "benchmark.yaml"
    benchmark_file.write_text(benchmark_content)
    return benchmark_file


@pytest.mark.short
def test_create_module_with_stage_for_initial_stage(tmp_path, sample_benchmark_yaml):
    """Test creating a module for an initial stage (no inputs)."""
    module_dir = tmp_path / "test-data-module"

    runner = CliRunner()
    result = runner.invoke(
        create_module,
        [
            str(module_dir),
            "--benchmark",
            str(sample_benchmark_yaml),
            "--for-stage",
            "data",
            "--non-interactive",
            "--name",
            "test-data-module",
            "--author-name",
            "Test Author",
            "--author-email",
            "test@example.com",
            "--entrypoint",
            "run.py",
        ],
    )

    assert result.exit_code == 0, f"Failed with: {result.output}"

    # Check that the module was created
    assert module_dir.exists()
    assert (module_dir / "run.py").exists()
    assert (module_dir / "omnibenchmark.yaml").exists()

    # Check that entrypoint has correct CLI parsing (no stage inputs for initial stage)
    entrypoint_content = (module_dir / "run.py").read_text()
    assert "--output_dir" in entrypoint_content
    assert "--name" in entrypoint_content
    # Should not have stage-specific inputs since it's an initial stage
    assert "--data.raw" not in entrypoint_content


@pytest.mark.short
def test_create_module_with_stage_for_methods_stage(tmp_path, sample_benchmark_yaml):
    """Test creating a module for methods stage with inputs."""
    module_dir = tmp_path / "test-method-module"

    runner = CliRunner()
    result = runner.invoke(
        create_module,
        [
            str(module_dir),
            "--benchmark",
            str(sample_benchmark_yaml),
            "--for-stage",
            "methods",
            "--non-interactive",
            "--name",
            "test-method-module",
            "--author-name",
            "Test Author",
            "--author-email",
            "test@example.com",
            "--entrypoint",
            "run.py",
        ],
    )

    assert result.exit_code == 0, f"Failed with: {result.output}"

    # Check that the module was created
    assert module_dir.exists()
    assert (module_dir / "run.py").exists()

    # Check that entrypoint has correct CLI parsing with stage inputs
    entrypoint_content = (module_dir / "run.py").read_text()
    assert "--output_dir" in entrypoint_content
    assert "--name" in entrypoint_content
    assert "--data.raw" in entrypoint_content
    assert "--data.meta" in entrypoint_content

    # Check that inputs use nargs='+' for multiple values
    assert "nargs='+'" in entrypoint_content

    # Check dest parameter converts dots to underscores
    assert "dest='data_raw'" in entrypoint_content
    assert "dest='data_meta'" in entrypoint_content


@pytest.mark.short
def test_create_module_with_stage_r_entrypoint(tmp_path, sample_benchmark_yaml):
    """Test creating an R module with stage inputs."""
    module_dir = tmp_path / "test-method-r-module"

    runner = CliRunner()
    result = runner.invoke(
        create_module,
        [
            str(module_dir),
            "--benchmark",
            str(sample_benchmark_yaml),
            "--for-stage",
            "metrics",
            "--non-interactive",
            "--name",
            "test-method-r-module",
            "--author-name",
            "Test Author",
            "--author-email",
            "test@example.com",
            "--entrypoint",
            "run.R",
        ],
    )

    assert result.exit_code == 0, f"Failed with: {result.output}"

    # Check that the module was created
    assert module_dir.exists()
    assert (module_dir / "run.R").exists()

    # Check that entrypoint has correct CLI parsing
    entrypoint_content = (module_dir / "run.R").read_text()
    assert "--output_dir" in entrypoint_content
    assert "--name" in entrypoint_content
    assert "--methods.result" in entrypoint_content
    assert "--data.raw" in entrypoint_content

    # Check R-specific syntax
    assert 'nargs="+"' in entrypoint_content
    assert 'dest="methods_result"' in entrypoint_content
    assert 'dest="data_raw"' in entrypoint_content


@pytest.mark.short
def test_create_module_with_stage_sh_entrypoint(tmp_path, sample_benchmark_yaml):
    """Test creating a shell script module with stage inputs."""
    module_dir = tmp_path / "test-method-sh-module"

    runner = CliRunner()
    result = runner.invoke(
        create_module,
        [
            str(module_dir),
            "--benchmark",
            str(sample_benchmark_yaml),
            "--for-stage",
            "methods",
            "--non-interactive",
            "--name",
            "test-method-sh-module",
            "--author-name",
            "Test Author",
            "--author-email",
            "test@example.com",
            "--entrypoint",
            "run.sh",
        ],
    )

    assert result.exit_code == 0, f"Failed with: {result.output}"

    # Check that the module was created
    assert module_dir.exists()
    assert (module_dir / "run.sh").exists()

    # Check that entrypoint has correct CLI parsing
    entrypoint_content = (module_dir / "run.sh").read_text()
    assert "--output_dir" in entrypoint_content
    assert "--name" in entrypoint_content
    assert "--data.raw" in entrypoint_content
    assert "--data.meta" in entrypoint_content

    # Check bash-specific array handling
    assert "DATA_RAW=()" in entrypoint_content
    assert "DATA_META=()" in entrypoint_content


@pytest.mark.short
def test_create_module_for_stage_requires_benchmark(tmp_path):
    """Test that --for-stage requires --benchmark."""
    module_dir = tmp_path / "test-module"

    runner = CliRunner()
    result = runner.invoke(
        create_module,
        [
            str(module_dir),
            "--for-stage",
            "methods",
            "--non-interactive",
            "--name",
            "test-module",
            "--author-name",
            "Test Author",
            "--author-email",
            "test@example.com",
        ],
    )

    assert result.exit_code != 0
    assert "--for-stage requires --benchmark" in result.output


@pytest.mark.short
def test_create_module_benchmark_requires_for_stage(tmp_path, sample_benchmark_yaml):
    """Test that --benchmark requires --for-stage."""
    module_dir = tmp_path / "test-module"

    runner = CliRunner()
    result = runner.invoke(
        create_module,
        [
            str(module_dir),
            "--benchmark",
            str(sample_benchmark_yaml),
            "--non-interactive",
            "--name",
            "test-module",
            "--author-name",
            "Test Author",
            "--author-email",
            "test@example.com",
        ],
    )

    assert result.exit_code != 0
    assert "--benchmark requires --for-stage" in result.output


@pytest.mark.short
def test_create_module_with_invalid_stage(tmp_path, sample_benchmark_yaml):
    """Test error handling for invalid stage ID."""
    module_dir = tmp_path / "test-module"

    runner = CliRunner()
    result = runner.invoke(
        create_module,
        [
            str(module_dir),
            "--benchmark",
            str(sample_benchmark_yaml),
            "--for-stage",
            "nonexistent",
            "--non-interactive",
            "--name",
            "test-module",
            "--author-name",
            "Test Author",
            "--author-email",
            "test@example.com",
        ],
    )

    assert result.exit_code != 0
    assert "Stage 'nonexistent' not found" in result.output
    assert "Available stages:" in result.output


@pytest.mark.short
def test_create_module_without_stage(tmp_path):
    """Test creating a module without stage-specific inputs (backward compatibility)."""
    module_dir = tmp_path / "test-generic-module"

    runner = CliRunner()
    result = runner.invoke(
        create_module,
        [
            str(module_dir),
            "--non-interactive",
            "--name",
            "test-generic-module",
            "--author-name",
            "Test Author",
            "--author-email",
            "test@example.com",
            "--entrypoint",
            "run.py",
        ],
    )

    assert result.exit_code == 0, f"Failed with: {result.output}"

    # Check that the module was created
    assert module_dir.exists()
    assert (module_dir / "run.py").exists()

    # Check that entrypoint has basic CLI parsing without stage inputs
    entrypoint_content = (module_dir / "run.py").read_text()
    assert "--output_dir" in entrypoint_content
    assert "--name" in entrypoint_content
    # Should have comment about adding custom inputs
    assert "Add your custom input arguments here" in entrypoint_content
