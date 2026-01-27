import pytest
from pathlib import Path

from tests.e2e.common import (
    run_standard_pipeline_test,
)


@pytest.fixture
def data_and_methods_config():
    """Get the path to the data and methods config."""
    config_path = Path(__file__).parent / "configs" / "01_data_and_methods.yaml"
    return config_path


@pytest.mark.e2e
def test_data_and_methods_pipeline(
    data_and_methods_config, tmp_path, bundled_repos, keep_files
):
    """Test a data and methods pipeline using the omnibenchmark CLI."""
    run_standard_pipeline_test(
        config_path=data_and_methods_config,
        config_filename="01_data_and_methods.yaml",
        test_name="01_data_and_methods",
        tmp_path=tmp_path,
        keep_files=keep_files,
        min_expected_files=4,  # We expect 4 *_data.json files (2 datasets + 2 method results)
        additional_cli_args=["-y"],
    )
