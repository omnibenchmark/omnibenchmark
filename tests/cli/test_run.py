import pytest

from click.testing import CliRunner
from pathlib import Path
from unittest.mock import patch, MagicMock
from omnibenchmark.cli.run import run_benchmark, run_module, validate_yaml
from omnibenchmark.model.validation import BenchmarkParseError

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
    extra arguments that are passed directly to nakemake.
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
            "--jobs",
            "6",
            "--cores",
            "24",
            "--default-resources",
            "mem_mb=4000",
            "runtime=600",
            "--verbose",
            "--printshellcmds",
            "--yes",
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


@pytest.mark.short
def test_run_benchmark_with_parse_error():
    """
    Test that run_benchmark handles BenchmarkParseError correctly
    and displays formatted error message.
    """
    benchmark_path = Path(data / "mock_benchmark.yaml").as_posix()

    parse_error = BenchmarkParseError("Invalid parameter type")
    parse_error.file_path = benchmark_path
    parse_error.line_number = 10

    with patch("omnibenchmark.cli.run.BenchmarkExecution") as mock_execution:
        mock_execution.side_effect = parse_error

        with patch("omnibenchmark.cli.run.pretty_print_parse_error") as mock_format:
            mock_format.return_value = "Formatted error message"

            runner = CliRunner()
            result = runner.invoke(
                run_benchmark,
                ["--benchmark", benchmark_path, "--local-storage", "--yes"],
            )

            # Verify the error was formatted
            mock_format.assert_called_once_with(parse_error)

            # Verify the command failed
            assert result.exit_code == 1


@pytest.mark.short
def test_run_benchmark_with_generic_exception():
    """
    Test that run_benchmark handles generic exceptions correctly.
    """
    benchmark_path = Path(data / "mock_benchmark.yaml").as_posix()

    with patch("omnibenchmark.cli.run.BenchmarkExecution") as mock_execution:
        mock_execution.side_effect = RuntimeError("Something went wrong")

        runner = CliRunner()
        result = runner.invoke(
            run_benchmark,
            ["--benchmark", benchmark_path, "--local-storage", "--yes"],
        )

        # Verify the command failed
        assert result.exit_code == 1


@pytest.mark.short
def test_run_module_with_parse_error():
    """
    Test that run_module handles BenchmarkParseError correctly.
    """
    benchmark_path = Path(data / "mock_benchmark.yaml").as_posix()

    parse_error = BenchmarkParseError("Invalid module configuration")
    parse_error.stage_id = "clustering"
    parse_error.module_id = "test_module"

    with patch("omnibenchmark.cli.run.BenchmarkExecution") as mock_execution:
        mock_execution.side_effect = parse_error

        with patch("omnibenchmark.cli.run.pretty_print_parse_error") as mock_format:
            mock_format.return_value = "Formatted module error"

            runner = CliRunner()
            result = runner.invoke(
                run_module,
                ["--benchmark", benchmark_path, "--module", "test_module"],
            )

            # Verify the error was formatted
            mock_format.assert_called_once_with(parse_error)

            # Verify the command failed
            assert result.exit_code == 1


@pytest.mark.short
def test_run_module_with_generic_exception():
    """
    Test that run_module handles generic exceptions correctly.
    """
    benchmark_path = Path(data / "mock_benchmark.yaml").as_posix()

    with patch("omnibenchmark.cli.run.BenchmarkExecution") as mock_execution:
        mock_execution.side_effect = ValueError("Invalid configuration")

        runner = CliRunner()
        result = runner.invoke(
            run_module,
            ["--benchmark", benchmark_path, "--module", "test_module"],
        )

        # Verify the command failed
        assert result.exit_code == 1


@pytest.mark.short
def test_validate_yaml_with_parse_error():
    """
    Test that validate_yaml handles BenchmarkParseError correctly.
    """
    benchmark_path = Path(data / "mock_benchmark.yaml").as_posix()

    parse_error = BenchmarkParseError("Validation failed")
    parse_error.file_path = benchmark_path

    with patch("omnibenchmark.cli.run.BenchmarkExecution") as mock_execution:
        mock_execution.side_effect = parse_error

        with patch("omnibenchmark.cli.run.pretty_print_parse_error") as mock_format:
            mock_format.return_value = "Formatted validation error"

            runner = CliRunner()
            result = runner.invoke(
                validate_yaml,
                ["--benchmark", benchmark_path],
            )

            # Verify the error was formatted
            mock_format.assert_called_once_with(parse_error)

            # Verify the command failed
            assert result.exit_code == 1


@pytest.mark.short
def test_validate_yaml_with_generic_exception():
    """
    Test that validate_yaml handles generic exceptions correctly.
    """
    benchmark_path = Path(data / "mock_benchmark.yaml").as_posix()

    with patch("omnibenchmark.cli.run.BenchmarkExecution") as mock_execution:
        mock_execution.side_effect = IOError("Cannot read file")

        runner = CliRunner()
        result = runner.invoke(
            validate_yaml,
            ["--benchmark", benchmark_path],
        )

        # Verify the command failed
        assert result.exit_code == 1


@pytest.mark.short
def test_validate_yaml_success():
    """
    Test that validate_yaml succeeds with a valid benchmark.
    """
    benchmark_path = Path(data / "mock_benchmark.yaml").as_posix()

    with patch("omnibenchmark.cli.run.BenchmarkExecution") as mock_execution:
        mock_execution.return_value = MagicMock()

        runner = CliRunner()
        result = runner.invoke(
            validate_yaml,
            ["--benchmark", benchmark_path],
        )

        # Verify the command succeeded
        assert result.exit_code == 0


@pytest.mark.short
def test_run_benchmark_with_top_level_field_parse_error():
    """
    Test that run_benchmark handles BenchmarkParseError for top-level fields
    (e.g., 'version', 'storage') and provides line context.
    """
    benchmark_path = Path(data / "mock_benchmark.yaml").as_posix()

    # Create a parse error for a top-level field like 'version'
    parse_error = BenchmarkParseError("Input should be a valid string")
    parse_error.file_path = benchmark_path
    parse_error.line_number = 3

    with patch("omnibenchmark.cli.run.BenchmarkExecution") as mock_execution:
        mock_execution.side_effect = parse_error

        with patch("omnibenchmark.cli.run.pretty_print_parse_error") as mock_format:
            mock_format.return_value = "Formatted top-level error"

            runner = CliRunner()
            result = runner.invoke(
                run_benchmark,
                ["--benchmark", benchmark_path, "--local-storage", "--yes"],
            )

            # Verify the error was formatted
            mock_format.assert_called_once_with(parse_error)

            # Verify the command failed
            assert result.exit_code == 1
