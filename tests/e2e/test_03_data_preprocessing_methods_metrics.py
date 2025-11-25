import pytest
from pathlib import Path


@pytest.fixture
def data_preprocessing_methods_metrics_config():
    """Get the path to the data, preprocessing, methods and metrics config."""
    config_path = (
        Path(__file__).parent / "configs" / "03_data_preprocessing_methods_metrics.yaml"
    )
    return config_path


@pytest.mark.e2e
def test_data_preprocessing_methods_metrics_pipeline(
    data_preprocessing_methods_metrics_config, tmp_path, bundled_repos, keep_files
):
    """Test a data, preprocessing, methods and metrics pipeline using the omnibenchmark CLI."""
    from tests.e2e.common import E2ETestRunner

    runner = E2ETestRunner(tmp_path, keep_files)
    config_file_in_tmp = runner.setup_test_environment(
        data_preprocessing_methods_metrics_config,
        "03_data_preprocessing_methods_metrics.yaml",
    )

    # Execute CLI with default args
    runner.execute_cli_command(config_file_in_tmp, ["--continue-on-error", "-y"])

    # Validate results
    runner.validate_results("03_data_preprocessing_methods_metrics")

    # Verify different file types separately
    runner.verify_output_file_count(2, "*_data.json")  # 2 data files
    runner.verify_output_file_count(4, "*_preprocessed*.json")  # 4 preprocessing files
    runner.verify_output_file_count(8, "*_method*.json")  # 8 method result files
