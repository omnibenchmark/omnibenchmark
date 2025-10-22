import pytest
from pathlib import Path

from tests.e2e.common import (
    E2ETestRunner,
    run_standard_pipeline_test,
    compare_pipeline_runs,
    extract_test_name_from_config,
)


@pytest.fixture
def data_modules_config():
    """Get the path to the data modules config."""
    config_path = Path(__file__).parent / "configs" / "00_data_modules.yaml"
    return config_path


@pytest.mark.e2e
def test_data_modules_pipeline(
    data_modules_config, tmp_path, bundled_repos, keep_files
):
    """Test a data modules pipeline using the omnibenchmark CLI.

    This is a black-box test - we execute the CLI and validate outputs.
    """
    run_standard_pipeline_test(
        config_path=data_modules_config,
        config_filename="00_data_modules.yaml",
        test_name="00_data_modules",
        tmp_path=tmp_path,
        keep_files=keep_files,
        min_expected_files=2,  # We expect 2 data files (D1_data.json, D2_data.json, excluding symlinks)
        additional_cli_args=["-y"],
    )


@pytest.mark.e2e
def test_data_modules_output_structure(
    data_modules_config, tmp_path, bundled_repos, keep_files
):
    """Test that the data modules pipeline creates the expected output directory structure."""
    run_standard_pipeline_test(
        config_path=data_modules_config,
        config_filename="00_data_modules.yaml",
        test_name="00_data_modules",
        tmp_path=tmp_path,
        keep_files=keep_files,
        min_expected_files=2,
        additional_cli_args=["-y"],
    )


@pytest.mark.e2e
def test_data_modules_idempotent(
    data_modules_config, tmp_path, bundled_repos, keep_files
):
    """Test that running the same workflow twice produces identical results."""
    # Use same tmp_path but different output directories
    first_runner = E2ETestRunner(tmp_path, keep_files)
    first_runner.out_dir = tmp_path / "out_first"

    second_runner = E2ETestRunner(tmp_path, keep_files)
    second_runner.out_dir = tmp_path / "out_second"

    config_filename = "00_data_modules.yaml"
    test_name = extract_test_name_from_config(config_filename)

    # Setup config file once
    config_file = first_runner.setup_test_environment(
        data_modules_config, config_filename
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


@pytest.mark.e2e
def test_data_modules_cli_validation(
    data_modules_config, tmp_path, bundled_repos, keep_files
):
    """Test that the CLI handles the benchmark config correctly."""
    # This test is functionally identical to the basic pipeline test
    # The CLI validation is inherent in the execution and result validation
    run_standard_pipeline_test(
        config_path=data_modules_config,
        config_filename="00_data_modules.yaml",
        test_name="00_data_modules",
        tmp_path=tmp_path,
        keep_files=keep_files,
        min_expected_files=2,
        additional_cli_args=["-y"],
    )
