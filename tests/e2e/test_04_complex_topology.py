import pytest
from pathlib import Path

from tests.e2e.common import (
    E2ETestRunner,
)


@pytest.fixture
def complex_topology_config():
    """Get the path to the complex topology config."""
    config_path = Path(__file__).parent / "configs" / "04_complex_topology.yaml"
    return config_path


@pytest.mark.e2e
def test_complex_topology_pipeline(
    complex_topology_config, tmp_path, bundled_repos, keep_files
):
    """Test a complex topology pipeline with multiple method groups using the omnibenchmark CLI."""
    runner = E2ETestRunner(tmp_path, keep_files)
    config_file_in_tmp = runner.setup_test_environment(
        complex_topology_config,
        "04_complex_topology.yaml",
    )

    # Execute CLI with default args
    runner.execute_cli_command(config_file_in_tmp, ["--continue-on-error", "-y"])

    # Validate results
    runner.validate_results("04_complex_topology")

    # Verify different file types for the complex topology:
    # - 4 data files (D1, D2, D3, D4)
    # - 4 preprocessing files (all datasets get normalized)
    # - Group 1: 4 results (2 primary datasets * 2 methods working on preprocessed)
    # - Group 1S: 4 results (2 secondary datasets * 2 methods working on preprocessed)
    # - Group 2: 2 results (2 primary datasets * 1 method working on raw)
    # - Group 3: 2 results (2 secondary datasets * 1 method working on raw)
    # Total method results: 4 + 4 + 2 + 2 = 12 results
    runner.verify_output_file_count(4, "*_data.json")  # 4 data files
    runner.verify_output_file_count(4, "*_normalized.json")  # 4 preprocessing files
    runner.verify_output_file_count(12, "*_method.json")  # 12 total method results


@pytest.mark.e2e
def test_complex_topology_idempotent(
    complex_topology_config, tmp_path, bundled_repos, keep_files
):
    """Test that running the complex topology workflow twice produces identical results."""
    from tests.e2e.common import (
        compare_pipeline_runs,
        extract_test_name_from_config,
    )

    # Use same tmp_path but different output directories
    first_runner = E2ETestRunner(tmp_path, keep_files)
    first_runner.out_dir = tmp_path / "out_first"

    second_runner = E2ETestRunner(tmp_path, keep_files)
    second_runner.out_dir = tmp_path / "out_second"

    config_filename = "04_complex_topology.yaml"
    test_name = extract_test_name_from_config(config_filename)

    # Setup config file once
    config_file = first_runner.setup_test_environment(
        complex_topology_config, config_filename
    )

    # First run
    first_runner.execute_cli_command(
        config_file, ["--continue-on-error", "-y"], debug_label="first run"
    )
    first_runner.validate_results(test_name)

    # Second run (same config, different output dir)
    second_runner.execute_cli_command(
        config_file, ["--continue-on-error", "-y"], debug_label="second run"
    )
    second_runner.validate_results(test_name)

    # Compare output directories
    compare_pipeline_runs(first_runner.out_dir, second_runner.out_dir)
