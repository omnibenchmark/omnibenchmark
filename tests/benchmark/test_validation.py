import tempfile
from pathlib import Path

import pytest

from omnibenchmark.benchmark.converter import LinkMLConverter
from omnibenchmark.benchmark.validation.validator import Validator


@pytest.mark.short
def test_validate_output_patterns_no_ambiguity():
    """Test that validate_output_patterns returns no errors when outputs are unique."""
    benchmark_file = "../data/Benchmark_001.yaml"
    benchmark_file_path = Path(__file__).parent / benchmark_file

    converter = LinkMLConverter(benchmark_file_path)
    errors = Validator.validate_output_patterns(converter)

    assert len(errors) == 0, f"Expected no errors but got: {errors}"


@pytest.mark.short
def test_validate_output_patterns_with_ambiguity():
    """Test that validate_output_patterns detects ambiguous output patterns."""
    # Create a temporary YAML file with ambiguous patterns
    ambiguous_yaml = """id: test_ambiguous
description: Test benchmark with ambiguous output patterns
version: 1.0
benchmarker: "Test User"
benchmark_yaml_spec: 0.01
software_backend: host
software_environments:
  env1:
    description: "test environment"
stages:
  - id: stage1
    modules:
      - id: module1
        name: "Module 1"
        software_environment: "env1"
        repository:
          url: https://github.com/test/repo.git
          commit: abc123
    outputs:
      - id: stage1.output
        path: "{dataset}.txt"

  - id: stage2
    modules:
      - id: module2
        name: "Module 2"
        software_environment: "env1"
        repository:
          url: https://github.com/test/repo.git
          commit: abc123
    outputs:
      - id: stage2.output
        path: "{dataset}.txt"
"""

    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
        f.write(ambiguous_yaml)
        temp_file = f.name

    try:
        converter = LinkMLConverter(Path(temp_file))
        errors = Validator.validate_output_patterns(converter)

        # Both stages output {dataset}.txt
        assert len(errors) > 0, "Expected to find ambiguous output patterns"

        # Check that the error message mentions the ambiguous pattern
        error_message = str(errors[0])
        assert (
            "{dataset}.txt" in error_message
        ), "Error should mention the ambiguous pattern"
        assert "stage1" in error_message, "Error should mention stage1"
        assert "stage2" in error_message, "Error should mention stage2"
        assert (
            "Ambiguous output pattern" in error_message
        ), "Error should mention ambiguity"
    finally:
        Path(temp_file).unlink()


@pytest.mark.short
def test_validate_output_patterns_different_directories():
    """Test that outputs with same filename but different directories are still flagged."""
    # Create a temporary YAML with same filename in different directories
    yaml_content = """id: test_dirs
description: Test benchmark with outputs in different directories
version: 1.0
benchmarker: "Test User"
benchmark_yaml_spec: 0.01
software_backend: host
software_environments:
  env1:
    description: "test environment"
stages:
  - id: stage1
    modules:
      - id: module1
        name: "Module 1"
        software_environment: "env1"
        repository:
          url: https://github.com/test/repo.git
          commit: abc123
    outputs:
      - id: stage1.output
        path: "dir1/{dataset}.txt"

  - id: stage2
    modules:
      - id: module2
        name: "Module 2"
        software_environment: "env1"
        repository:
          url: https://github.com/test/repo.git
          commit: abc123
    outputs:
      - id: stage2.output
        path: "dir2/{dataset}.txt"
"""

    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
        f.write(yaml_content)
        temp_file = f.name

    try:
        converter = LinkMLConverter(Path(temp_file))
        errors = Validator.validate_output_patterns(converter)

        # Note: The current implementation only checks base filenames
        # If two stages have outputs like "dir1/{dataset}.txt" and "dir2/{dataset}.txt",
        # they will be flagged as ambiguous because os.path.basename returns the same pattern
        # This is actually correct behavior for Snakemake, as it can create ambiguous rules
        assert (
            len(errors) > 0
        ), "Expected ambiguous patterns even in different directories"
    finally:
        Path(temp_file).unlink()


@pytest.mark.short
def test_validate_output_patterns_multiple_ambiguities():
    """Test detection of multiple different ambiguous patterns."""
    yaml_content = """id: test_multiple
description: Test benchmark with multiple ambiguous patterns
version: 1.0
benchmarker: "Test User"
benchmark_yaml_spec: 0.01
software_backend: host
software_environments:
  env1:
    description: "test environment"
stages:
  - id: stage1
    modules:
      - id: module1
        name: "Module 1"
        software_environment: "env1"
        repository:
          url: https://github.com/test/repo.git
          commit: abc123
    outputs:
      - id: stage1.output1
        path: "{dataset}.txt"
      - id: stage1.output2
        path: "{dataset}.json"

  - id: stage2
    modules:
      - id: module2
        name: "Module 2"
        software_environment: "env1"
        repository:
          url: https://github.com/test/repo.git
          commit: abc123
    outputs:
      - id: stage2.output1
        path: "{dataset}.txt"
      - id: stage2.output2
        path: "{dataset}.json"
"""

    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
        f.write(yaml_content)
        temp_file = f.name

    try:
        converter = LinkMLConverter(Path(temp_file))
        errors = Validator.validate_output_patterns(converter)

        # Should find two ambiguous patterns: {dataset}.txt and {dataset}.json
        assert len(errors) == 2, f"Expected 2 ambiguous patterns but got {len(errors)}"
    finally:
        Path(temp_file).unlink()


@pytest.mark.short
def test_validate_output_patterns_single_stage_multiple_outputs():
    """Test that multiple outputs within the same stage don't trigger ambiguity errors."""
    # A single stage with multiple outputs should not be flagged as ambiguous
    benchmark_file = "../data/Benchmark_001.yaml"
    benchmark_file_path = Path(__file__).parent / benchmark_file

    converter = LinkMLConverter(benchmark_file_path)

    # Check that at least one stage has multiple outputs
    stages = converter.get_stages()
    has_stage_with_multiple_outputs = False
    for stage in stages.values():
        if stage.outputs and len(stage.outputs) > 1:
            has_stage_with_multiple_outputs = True
            break

    # Benchmark_001 has stages with multiple outputs (data stage has 3 outputs)
    assert (
        has_stage_with_multiple_outputs
    ), "Test file should have a stage with multiple outputs"

    errors = Validator.validate_output_patterns(converter)

    # Should have no errors for unique patterns across stages
    assert len(errors) == 0


@pytest.mark.short
def test_validate_output_patterns_error_message_format():
    """Test that error messages contain helpful information for fixing the issue."""
    yaml_content = """id: test_error_format
description: Test benchmark for error message format
version: 1.0
benchmarker: "Test User"
benchmark_yaml_spec: 0.01
software_backend: host
software_environments:
  env1:
    description: "test environment"
stages:
  - id: rawdata
    modules:
      - id: module1
        name: "Module 1"
        software_environment: "env1"
        repository:
          url: https://github.com/test/repo.git
          commit: abc123
    outputs:
      - id: rawdata.output
        path: "{dataset}.ad"

  - id: simulate
    modules:
      - id: module2
        name: "Module 2"
        software_environment: "env1"
        repository:
          url: https://github.com/test/repo.git
          commit: abc123
    outputs:
      - id: simulate.output
        path: "{dataset}.ad"
"""

    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
        f.write(yaml_content)
        temp_file = f.name

    try:
        converter = LinkMLConverter(Path(temp_file))
        errors = Validator.validate_output_patterns(converter)

        assert len(errors) > 0, "Expected to find ambiguous patterns"

        error_message = str(errors[0])

        # Check that the error message is helpful
        assert "stage" in error_message.lower(), "Error should mention stages"
        assert "output" in error_message.lower(), "Error should mention outputs"

        # Check that it provides actionable suggestions
        assert (
            "suggestion" in error_message.lower() or "modify" in error_message.lower()
        ), "Error should provide suggestions for fixing the issue"
    finally:
        Path(temp_file).unlink()
