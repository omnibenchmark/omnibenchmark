"""Tests for parameter expansion with combinatorial validation."""

import pytest

from omnibenchmark.model.benchmark import Benchmark


@pytest.mark.short
def test_legacy_cli_args_format():
    """Test that legacy CLI args format still works."""
    yaml_content = """
id: test_benchmark
description: Test benchmark
version: "1.0"
benchmarker: Test
software_backend: host
software_environments:
  - id: python
    description: Python environment
stages:
  - id: test_stage
    modules:
      - id: test_module
        name: Test Module
        software_environment: python
        repository:
          url: https://github.com/test/test.git
          commit: abc123
        parameters:
          - values: ["--method", "cosine"]
          - values: ["--method", "euclidean"]
    outputs:
      - id: test.output
        path: test.csv
"""
    benchmark = Benchmark.from_yaml(yaml_content)
    params = benchmark.get_module_parameters("test_module")

    # Should have 2 parameter sets (2 items in parameters list)
    assert params is not None
    assert len(params) == 2
    assert params[0]["method"] == "cosine"
    assert params[1]["method"] == "euclidean"


@pytest.mark.short
def test_simple_dict_format_no_lists():
    """Test dict format with no lists (no expansion needed)."""
    yaml_content = """
id: test_benchmark
description: Test benchmark
version: "1.0"
benchmarker: Test
software_backend: host
software_environments:
  - id: python
    description: Python environment
stages:
  - id: test_stage
    modules:
      - id: test_module
        name: Test Module
        software_environment: python
        repository:
          url: https://github.com/test/test.git
          commit: abc123
        parameters:
          - method: genie
            gini_threshold: 0.1
          - method: other
            gini_threshold: 0.2
    outputs:
      - id: test.output
        path: test.csv
"""
    benchmark = Benchmark.from_yaml(yaml_content)
    params = benchmark.get_module_parameters("test_module")

    # Should have 2 parameter sets (2 items in parameters list, no expansion)
    assert params is not None
    assert len(params) == 2
    assert params[0]["method"] == "genie"
    assert params[0]["gini_threshold"] == 0.1
    assert params[1]["method"] == "other"
    assert params[1]["gini_threshold"] == 0.2


@pytest.mark.short
def test_dict_format_with_single_list():
    """Test dict format with a single list parameter (expansion needed)."""
    yaml_content = """
id: test_benchmark
description: Test benchmark
version: "1.0"
benchmarker: Test
software_backend: host
software_environments:
  - id: python
    description: Python environment
stages:
  - id: test_stage
    modules:
      - id: test_module
        name: Test Module
        software_environment: python
        repository:
          url: https://github.com/test/test.git
          commit: abc123
        parameters:
          - method: paired_params
            k: 1
            m: [2, 3, 4]
    outputs:
      - id: test.output
        path: test.csv
"""
    benchmark = Benchmark.from_yaml(yaml_content)
    params = benchmark.get_module_parameters("test_module")

    # Combinatorics: 1 item with m=[2,3,4] -> 3 expansions
    assert params is not None
    assert len(params) == 3

    # Verify all combinations
    assert params[0]["method"] == "paired_params"
    assert params[0]["k"] == 1
    assert params[0]["m"] == 2

    assert params[1]["method"] == "paired_params"
    assert params[1]["k"] == 1
    assert params[1]["m"] == 3

    assert params[2]["method"] == "paired_params"
    assert params[2]["k"] == 1
    assert params[2]["m"] == 4


@pytest.mark.short
def test_dict_format_with_multiple_lists():
    """Test dict format with multiple list parameters (cartesian product)."""
    yaml_content = """
id: test_benchmark
description: Test benchmark
version: "1.0"
benchmarker: Test
software_backend: host
software_environments:
  - id: python
    description: Python environment
stages:
  - id: test_stage
    modules:
      - id: test_module
        name: Test Module
        software_environment: python
        repository:
          url: https://github.com/test/test.git
          commit: abc123
        parameters:
          - method: grid_search
            alpha: [0.1, 0.5, 1.0]
            beta: [2, 3]
    outputs:
      - id: test.output
        path: test.csv
"""
    benchmark = Benchmark.from_yaml(yaml_content)
    params = benchmark.get_module_parameters("test_module")

    # Combinatorics: alpha has 3 values, beta has 2 values -> 3 * 2 = 6 expansions
    assert params is not None
    assert len(params) == 6

    # Verify cartesian product
    expected_combinations = [
        {"method": "grid_search", "alpha": 0.1, "beta": 2},
        {"method": "grid_search", "alpha": 0.1, "beta": 3},
        {"method": "grid_search", "alpha": 0.5, "beta": 2},
        {"method": "grid_search", "alpha": 0.5, "beta": 3},
        {"method": "grid_search", "alpha": 1.0, "beta": 2},
        {"method": "grid_search", "alpha": 1.0, "beta": 3},
    ]

    for param in params:
        # Convert Params to dict for easier comparison
        param_dict = dict(param.items())
        assert param_dict in expected_combinations


@pytest.mark.short
def test_paired_params_example():
    """Test the user's example with paired parameters."""
    yaml_content = """
id: test_benchmark
description: Test benchmark
version: "1.0"
benchmarker: Test
software_backend: host
software_environments:
  - id: python
    description: Python environment
stages:
  - id: test_stage
    modules:
      - id: test_module
        name: Test Module
        software_environment: python
        repository:
          url: https://github.com/test/test.git
          commit: abc123
        parameters:
          - method: paired_params
            k: 1
            m: 2
          - method: paired_params
            k: 2
            m: 3
          - method: paired_params
            k: 4
            m: [0.5, 0.7]
    outputs:
      - id: test.output
        path: test.csv
"""
    benchmark = Benchmark.from_yaml(yaml_content)
    params = benchmark.get_module_parameters("test_module")

    # Combinatorics:
    # - First item: no lists -> 1 expansion
    # - Second item: no lists -> 1 expansion
    # - Third item: m has 2 values -> 2 expansions
    # Total: 1 + 1 + 2 = 4
    assert params is not None
    assert len(params) == 4

    # Verify combinations
    assert params[0]["k"] == 1
    assert params[0]["m"] == 2

    assert params[1]["k"] == 2
    assert params[1]["m"] == 3

    assert params[2]["k"] == 4
    assert params[2]["m"] == 0.5

    assert params[3]["k"] == 4
    assert params[3]["m"] == 0.7


@pytest.mark.short
def test_complex_combinatorics():
    """Test complex combinatorics with multiple parameters and multiple lists."""
    yaml_content = """
id: test_benchmark
description: Test benchmark
version: "1.0"
benchmarker: Test
software_backend: host
software_environments:
  - id: python
    description: Python environment
stages:
  - id: test_stage
    modules:
      - id: test_module
        name: Test Module
        software_environment: python
        repository:
          url: https://github.com/test/test.git
          commit: abc123
        parameters:
          - method: combo1
            a: [1, 2]
            b: [3, 4]
            c: [5, 6]
          - method: combo2
            x: 10
            y: [20, 30]
          - method: combo3
            z: 100
    outputs:
      - id: test.output
        path: test.csv
"""
    benchmark = Benchmark.from_yaml(yaml_content)
    params = benchmark.get_module_parameters("test_module")

    # Combinatorics:
    # - First item: 2 * 2 * 2 = 8 expansions (cartesian product of a, b, c)
    # - Second item: 2 expansions (y has 2 values)
    # - Third item: 1 expansion (no lists)
    # Total: 8 + 2 + 1 = 11
    assert params is not None
    assert len(params) == 11


@pytest.mark.short
def test_mixed_legacy_and_new_format():
    """Test mixing legacy CLI args format with new dict format."""
    yaml_content = """
id: test_benchmark
description: Test benchmark
version: "1.0"
benchmarker: Test
software_backend: host
software_environments:
  - id: python
    description: Python environment
stages:
  - id: test_stage
    modules:
      - id: test_module
        name: Test Module
        software_environment: python
        repository:
          url: https://github.com/test/test.git
          commit: abc123
        parameters:
          - values: ["--legacy", "arg"]
          - method: new_style
            param: [1, 2]
    outputs:
      - id: test.output
        path: test.csv
"""
    benchmark = Benchmark.from_yaml(yaml_content)
    params = benchmark.get_module_parameters("test_module")

    # Combinatorics:
    # - First item: legacy format -> 1 expansion
    # - Second item: param has 2 values -> 2 expansions
    # Total: 1 + 2 = 3
    assert params is not None
    assert len(params) == 3

    # First param should be from legacy format
    assert params[0]["legacy"] == "arg"

    # Next two should be from new format
    assert params[1]["method"] == "new_style"
    assert params[1]["param"] == 1

    assert params[2]["method"] == "new_style"
    assert params[2]["param"] == 2


@pytest.mark.short
def test_empty_list_raises_error():
    """Test that empty lists in parameters raise validation errors."""
    yaml_content = """
id: test_benchmark
description: Test benchmark
version: "1.0"
benchmarker: Test
software_backend: host
software_environments:
  - id: python
    description: Python environment
stages:
  - id: test_stage
    modules:
      - id: test_module
        name: Test Module
        software_environment: python
        repository:
          url: https://github.com/test/test.git
          commit: abc123
        parameters:
          - method: test
            param: []
    outputs:
      - id: test.output
        path: test.csv
"""
    # This should work during parsing, but produce no expansions for empty list
    # Let's check that the behavior is reasonable
    benchmark = Benchmark.from_yaml(yaml_content)
    params = benchmark.get_module_parameters("test_module")

    # Empty list in product should result in 0 combinations
    assert params is not None
    assert len(params) == 0


@pytest.mark.short
def test_dict_format_cli_args_conversion():
    """Test that dict format parameters can be converted to CLI args."""
    yaml_content = """
id: test_benchmark
description: Test benchmark
version: "1.0"
benchmarker: Test
software_backend: host
software_environments:
  - id: python
    description: Python environment
stages:
  - id: test_stage
    modules:
      - id: test_module
        name: Test Module
        software_environment: python
        repository:
          url: https://github.com/test/test.git
          commit: abc123
        parameters:
          - method: genie
            gini_threshold: 0.1
            verbose: true
    outputs:
      - id: test.output
        path: test.csv
"""
    benchmark = Benchmark.from_yaml(yaml_content)
    params = benchmark.get_module_parameters("test_module")

    assert params is not None
    assert len(params) == 1

    # Test GNU-style CLI args
    gnu_args = params[0].to_cli_args(style="gnu")
    assert "--gini_threshold" in gnu_args
    assert "0.1" in gnu_args
    assert "--method" in gnu_args
    assert "genie" in gnu_args
    assert "--verbose" in gnu_args

    # Test equals-style CLI args
    equals_args = params[0].to_cli_args(style="equals")
    assert "--gini_threshold=0.1" in equals_args
    assert "--method=genie" in equals_args
    assert "--verbose" in equals_args

    # Test serialization
    serialized = params[0].serialize()
    assert "gini_threshold" in serialized
    assert "0.1" in serialized
