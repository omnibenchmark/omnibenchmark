import pytest

from click.testing import CliRunner
from pathlib import Path
from unittest.mock import patch, MagicMock
from omnibenchmark.cli.run import run_benchmark

data = Path(__file__).parent.parent / "data"


@pytest.fixture
def mock_validate_benchmark():
    with patch("omnibenchmark.cli.run.validate_benchmark") as mock:
        mock.return_value = MagicMock()
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
    mock_validate_benchmark, mock_workflow_run_workflow, mock_click_confirm
):
    """
    Test that if we call run benchmark with --update and a valid benchmark file
    we get the click interactive confirmation prompt.
    """
    runner = CliRunner()
    benchmark_path = Path(data / "mock_benchmark.yaml").as_posix()

    runner = CliRunner()
    result = runner.invoke(
        run_benchmark,
        ["--benchmark", benchmark_path, "--cores", "2", "--update", "--local"],
    )

    # Ensure click.confirm is NOT called
    mock_click_confirm.assert_called_once_with(
        "Are you sure you want to re-run the entire workflow?", abort=True
    )

    # Ensure workflow.run_workflow is called
    mock_workflow_run_workflow.assert_called_once()
    assert result.exit_code == 0


@pytest.mark.short
def test_run_benchmark_with_yes(
    mock_validate_benchmark, mock_workflow_run_workflow, mock_click_confirm
):
    """
    Test that if we call run benchmark with -k, a valid benchmark file and the --yes flag
    we do not get the click interactive confirmation prompt, and run_workflow is called.
    """
    benchmark_path = Path(data / "mock_benchmark.yaml").as_posix()

    runner = CliRunner()
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
