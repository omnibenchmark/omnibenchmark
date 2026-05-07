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
    runner.execute_cli_command(config_file_in_tmp, ["--continue-on-error"])

    # Validate results
    runner.validate_results("03_data_preprocessing_methods_metrics")

    # Verify different file types separately
    runner.verify_output_file_count(2, "*_data.json")  # 2 data files
    runner.verify_output_file_count(4, "*_preprocessed*.json")  # 4 preprocessing files
    runner.verify_output_file_count(8, "*_method*.json")  # 8 method result files


@pytest.mark.e2e
def test_until_preprocessing_truncates_pipeline(
    data_preprocessing_methods_metrics_config, tmp_path, bundled_repos, keep_files
):
    """`--until preprocessing` runs data + preprocessing only; methods/metrics pruned."""
    from tests.e2e.common import (
        E2ETestRunner,
        filter_files_excluding_symlinked_dirs,
    )

    runner = E2ETestRunner(tmp_path, keep_files)
    config_file_in_tmp = runner.setup_test_environment(
        data_preprocessing_methods_metrics_config,
        "03_data_preprocessing_methods_metrics.yaml",
    )

    runner.execute_cli_command(
        config_file_in_tmp,
        ["--continue-on-error", "--until", "preprocessing"],
    )

    # data + preprocessing outputs are produced
    data_files = filter_files_excluding_symlinked_dirs(runner.out_dir, "*_data.json")
    preproc_files = filter_files_excluding_symlinked_dirs(
        runner.out_dir, "*_preprocessed*.json"
    )
    assert len(data_files) >= 2, f"expected ≥2 data files, got {len(data_files)}"
    assert (
        len(preproc_files) >= 4
    ), f"expected ≥4 preprocessing files, got {len(preproc_files)}"

    # methods stage and metric collector outputs must NOT be produced
    method_files = filter_files_excluding_symlinked_dirs(
        runner.out_dir, "*_method*.json"
    )
    metrics_files = filter_files_excluding_symlinked_dirs(
        runner.out_dir, "metrics.json"
    )
    assert method_files == [], f"unexpected method outputs: {method_files}"
    assert metrics_files == [], f"unexpected metrics outputs: {metrics_files}"


@pytest.mark.e2e
def test_until_unknown_stage_errors(
    data_preprocessing_methods_metrics_config, tmp_path, bundled_repos
):
    """`--until <unknown>` exits non-zero with a helpful message."""
    from tests.cli.cli_setup import OmniCLISetup

    runner_tmp = tmp_path / "out"
    config_in_tmp = tmp_path / "03_data_preprocessing_methods_metrics.yaml"
    import shutil

    tmp_path.mkdir(parents=True, exist_ok=True)
    shutil.copy2(data_preprocessing_methods_metrics_config, config_in_tmp)

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                str(config_in_tmp),
                "--out-dir",
                str(runner_tmp),
                "--until",
                "no_such_stage",
            ],
            cwd=str(tmp_path),
        )

    assert result.returncode != 0, "expected non-zero exit"
    combined = (result.stdout or "") + (result.stderr or "")
    assert "no_such_stage" in combined
