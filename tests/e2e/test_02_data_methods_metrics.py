import pytest
from pathlib import Path

from tests.e2e.common import (
    run_standard_pipeline_test,
)


@pytest.fixture
def data_methods_metrics_config():
    """Get the path to the data methods and metrics config."""
    config_path = Path(__file__).parent / "configs" / "02_data_methods_metrics.yaml"
    return config_path


@pytest.mark.e2e
def test_data_methods_metrics_pipeline(
    data_methods_metrics_config, tmp_path, bundled_repos, keep_files
):
    """Test a data, methods and metrics pipeline using the omnibenchmark CLI."""
    run_standard_pipeline_test(
        config_path=data_methods_metrics_config,
        config_filename="02_data_methods_metrics.yaml",
        test_name="02_data_methods_metrics",
        tmp_path=tmp_path,
        keep_files=keep_files,
        min_expected_files=4,  # We expect 4 *_data.json files (2 datasets + 2 method results, metrics.json not counted)
        additional_cli_args=["-y"],
    )
