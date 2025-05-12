"""
Unit tests for omnibenchmark.cli.run.

At this level, we're interested in catching logical behaviour _before_ passing the
execution to the engine (snakemake.)

- They should all run fast.
- Do not actually call snakemake, use pytest mocks.

"""

import logging
import pytest

from click.testing import CliRunner
from pathlib import Path
from unittest.mock import patch, MagicMock

from omnibenchmark.benchmark import Benchmark
from omnibenchmark.cli.run import run_benchmark

data = Path(__file__).parent.parent / "data"


@pytest.fixture
def mock_validate_benchmark():
    with patch("omnibenchmark.cli.run.validate_benchmark") as mock:
        # Create a real Benchmark object with a mock YAML path
        mock_benchmark = MagicMock(spec=Benchmark)
        mock_benchmark.out_dir = None

        # Define the set_out_dir method
        def set_out_dir(path):
            mock_benchmark.out_dir = path

        # Attach the method to the mock
        mock_benchmark.set_out_dir = set_out_dir

        mock.return_value = mock_benchmark
        yield mock


@pytest.fixture
def mock_workflow_run_workflow():
    with patch("omnibenchmark.cli.run.SnakemakeEngine.run_workflow") as mock:
        mock.return_value = True
        yield mock


@pytest.fixture
def mock_click_confirm():
    with patch("click.confirm", return_value=True) as mock:
        yield mock


@pytest.mark.short
def test_run_benchmark_without_yes(
    mock_validate_benchmark, mock_workflow_run_workflow, mock_click_confirm, caplog
):
    """
    Test that if we call run benchmark with --update and a valid benchmark file
    we get the click interactive confirmation prompt.
    """
    benchmark_path = Path(data / "mock_benchmark.yaml").as_posix()

    runner = CliRunner()

    # Set log level to capture logs
    with caplog.at_level(logging.DEBUG):
        result = runner.invoke(
            run_benchmark,
            ["--benchmark", benchmark_path, "--cores", "2", "--update", "--local"],
        )

    # Ensure click.confirm is called
    mock_click_confirm.assert_called_once_with(
        "Are you sure you want to re-run the entire workflow?", abort=True
    )

    # Print captured logs for debugging
    assert "Running benchmark" in caplog.text

    # Ensure workflow.run_workflow is called
    mock_workflow_run_workflow.assert_called_once()
    assert result.exit_code == 0


@pytest.mark.short
def test_run_benchmark_with_yes(
    mock_validate_benchmark, mock_workflow_run_workflow, mock_click_confirm, caplog
):
    """
    Test that if we call run benchmark with -k, a valid benchmark file and the --yes flag
    we do not get the click interactive confirmation prompt, and run_workflow is called.
    """
    benchmark_path = Path(data / "mock_benchmark.yaml").as_posix()

    runner = CliRunner()
    # Set log level to capture logs
    with caplog.at_level(logging.DEBUG):
        result = runner.invoke(
            run_benchmark,
            [
                "--benchmark",
                benchmark_path,
                "-k",
                "--cores",
                "2",
                "--local",
                "--yes",
            ],
        )

    mock_click_confirm.assert_not_called()

    # Ensure workflow.run_workflow is called
    mock_workflow_run_workflow.assert_called_once()
    assert result.exit_code == 0


@pytest.mark.short
def test_run_benchmark_with_out_dir(
    mock_validate_benchmark, mock_workflow_run_workflow, caplog
):
    """
    Test that if we call run benchmark with --out-dir, the Benchmark object
    passed to run_workflow is passed the correct out_dir set and validate_benchmark
    is called with the correct out_dir.

    There are other side effects in this feature (setting of environment variable)
    that are too messy to test here, but please be aware of those if you work on
    these in the future.
    """
    benchmark_path = Path(data / "mock_benchmark.yaml").as_posix()
    out_dir = "custom_out"

    runner = CliRunner()
    # Set log level to capture logs
    with caplog.at_level(logging.INFO):
        result = runner.invoke(
            run_benchmark,
            [
                "--benchmark",
                benchmark_path,
                "--local",
                "--out-dir",
                out_dir,
            ],
        )

    # Ensure validate_benchmark is called with the correct out_dir
    mock_validate_benchmark.assert_called_once_with(benchmark_path, out_dir)

    # Ensure workflow.run_workflow is called
    mock_workflow_run_workflow.assert_called_once()

    # Extract the arguments passed to run_workflow
    workflow_args = mock_workflow_run_workflow.call_args[1]

    # Assert that the out_dir is correctly passed to run_workflow
    assert workflow_args["out_dir"] == out_dir
    assert result.exit_code == 0
