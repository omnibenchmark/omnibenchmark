import pytest

from click.testing import CliRunner
from pathlib import Path
from unittest.mock import patch, MagicMock
from omnibenchmark.cli.run import run_benchmark

data = Path(__file__).parent.parent / "data"


@pytest.fixture
def mock_benchmark_execution():
    with patch("omnibenchmark.cli.run.BenchmarkExecution") as mock:
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
    mock_benchmark_execution, mock_workflow_run_workflow, mock_click_confirm
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
        ["--benchmark", benchmark_path, "--cores", "2", "--update", "--local-storage"],
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
    mock_benchmark_execution, mock_workflow_run_workflow, mock_click_confirm
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
            "--local-storage",
            "--yes",
        ],
    )

    mock_click_confirm.assert_not_called()

    # Ensure workflow.run_workflow is called
    mock_workflow_run_workflow.assert_called_once()
    assert result.exit_code == 0


@pytest.mark.short
def test_run_benchmark_with_slurm_executor(
    mock_benchmark_execution, mock_workflow_run_workflow, mock_click_confirm
):
    """
    Test that the benchmark runs with the SLURM executor and
    extra arguments that are passed directly to Snakemake using the -- separator.
    """
    benchmark_path = Path(data / "mock_benchmark.yaml").as_posix()

    runner = CliRunner()
    result = runner.invoke(
        run_benchmark,
        [
            "--benchmark",
            benchmark_path,
            "-k",
            "--executor",
            "slurm",
            "--cores",
            "24",
            "--yes",
            "--out",
            "foo",
            "--",
            "--jobs",
            "6",
            "--default-resources",
            "mem_mb=4000",
            "runtime=600",
            "--verbose",
            "--printshellcmds",
        ],
    )

    mock_click_confirm.assert_not_called()

    # Ensure workflow.run_workflow is called
    mock_workflow_run_workflow.assert_called_once()
    args, kwargs = mock_workflow_run_workflow.call_args
    assert kwargs["executor"] == "slurm"
    assert kwargs["jobs"] == "6"
    assert kwargs["cores"] == 24
    assert kwargs["default-resources"] == ["mem_mb=4000", "runtime=600"]
    assert kwargs["continue_on_error"] is True
    assert kwargs["verbose"] is True
    assert kwargs["printshellcmds"] is True

    assert result.exit_code == 0
