import pytest

from click.testing import CliRunner
from pathlib import Path
from unittest.mock import patch, MagicMock
from omnibenchmark.cli.run import run
from omnibenchmark.cli.validate import validate_plan
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
def test_run_without_yes(
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
        run,
        [benchmark_path, "--cores", "2", "--update"],
    )

    # Ensure click.confirm is NOT called
    mock_click_confirm.assert_called_once_with(
        "Are you sure you want to re-run the entire workflow?", abort=True
    )

    # Ensure workflow.run_workflow is called
    mock_workflow_run_workflow.assert_called_once()
    assert result.exit_code == 0


@pytest.mark.short
def test_run_with_yes(
    mock_benchmark_execution, mock_workflow_run_workflow, mock_click_confirm
):
    """
    Test that if we call run benchmark with -k, a valid benchmark file and the --yes flag
    we do not get the click interactive confirmation prompt, and run_workflow is called.
    """
    benchmark_path = Path(data / "mock_benchmark.yaml").as_posix()

    runner = CliRunner()
    result = runner.invoke(
        run,
        [
            benchmark_path,
            "-k",
            "--cores",
            "2",
            "--yes",
        ],
    )

    mock_click_confirm.assert_not_called()

    # Ensure workflow.run_workflow is called
    mock_workflow_run_workflow.assert_called_once()
    assert result.exit_code == 0


@pytest.mark.short
def test_run_with_slurm_executor(
    mock_benchmark_execution, mock_workflow_run_workflow, mock_click_confirm
):
    """
    Test that the benchmark runs with the SLURM executor and
    extra arguments that are passed directly to Snakemake using the -- separator.
    """
    benchmark_path = Path(data / "mock_benchmark.yaml").as_posix()

    runner = CliRunner()
    result = runner.invoke(
        run,
        [
            benchmark_path,
            "-k",
            "--executor",
            "slurm",
            "--cores",
            "24",
            "--yes",
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


@pytest.mark.short
def test_run_with_parse_error():
    """
    Test that run handles BenchmarkParseError correctly
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
                run,
                [benchmark_path, "--yes"],
            )

            # Verify the error was formatted
            mock_format.assert_called_once_with(parse_error)

            # Verify the command failed
            assert result.exit_code == 1


@pytest.mark.short
def test_run_with_generic_exception():
    """
    Test that run handles generic exceptions correctly.
    """
    benchmark_path = Path(data / "mock_benchmark.yaml").as_posix()

    with patch("omnibenchmark.cli.run.BenchmarkExecution") as mock_execution:
        mock_execution.side_effect = RuntimeError("Something went wrong")

        runner = CliRunner()
        result = runner.invoke(
            run,
            [benchmark_path, "--yes"],
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
                run,
                [benchmark_path, "--module", "test_module"],
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
            run,
            [benchmark_path, "--module", "test_module"],
        )

        # Verify the command failed
        assert result.exit_code == 1


@pytest.mark.short
def test_validate_yaml_with_parse_error():
    """
    Test that validate plan handles BenchmarkParseError correctly.
    """
    benchmark_path = Path(data / "mock_benchmark.yaml").as_posix()

    parse_error = BenchmarkParseError("Validation failed")
    parse_error.file_path = benchmark_path

    with patch("omnibenchmark.cli.validate.BenchmarkModel") as mock_model:
        mock_model.from_yaml.side_effect = parse_error

        runner = CliRunner()
        result = runner.invoke(
            validate_plan,
            [benchmark_path],
        )

        # Verify the command failed
        assert result.exit_code == 1


@pytest.mark.short
def test_validate_yaml_with_generic_exception():
    """
    Test that validate plan handles generic exceptions correctly.
    """
    benchmark_path = Path(data / "mock_benchmark.yaml").as_posix()

    with patch("omnibenchmark.cli.validate.BenchmarkModel") as mock_model:
        mock_model.from_yaml.side_effect = IOError("Cannot read file")

        runner = CliRunner()
        result = runner.invoke(
            validate_plan,
            [benchmark_path],
        )

        # Verify the command failed
        assert result.exit_code == 1


@pytest.mark.short
def test_validate_yaml_success():
    """
    Test that validate plan succeeds with a valid benchmark.
    """
    benchmark_path = Path(data / "mock_benchmark.yaml").as_posix()

    runner = CliRunner()
    result = runner.invoke(
        validate_plan,
        [benchmark_path],
    )

    # Verify the command succeeded
    assert result.exit_code == 0


@pytest.mark.short
def test_run_with_top_level_field_parse_error():
    """
    Test that run handles BenchmarkParseError for top-level fields
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
                run,
                [benchmark_path, "--yes"],
            )

            # Verify the error was formatted
            mock_format.assert_called_once_with(parse_error)

            # Verify the command failed
            assert result.exit_code == 1


@pytest.mark.short
def test_run_continue_on_error_without_yes(
    mock_benchmark_execution, mock_workflow_run_workflow, mock_click_confirm
):
    """
    Test that --continue-on-error without --yes prompts for confirmation.
    """
    benchmark_path = Path(data / "mock_benchmark.yaml").as_posix()

    runner = CliRunner()
    result = runner.invoke(
        run,
        [benchmark_path, "--continue-on-error"],
    )

    # Ensure click.confirm was called for continue-on-error
    mock_click_confirm.assert_called_once_with(
        "Are you sure you want to run the full benchmark even if some jobs fail?",
        abort=True,
    )

    mock_workflow_run_workflow.assert_called_once()
    assert result.exit_code == 0


@pytest.mark.short
def test_run_with_remote_storage(mock_benchmark_execution, mock_workflow_run_workflow):
    """
    Test that --use-remote-storage flag sets up remote storage correctly.
    """
    benchmark_path = Path(data / "mock_benchmark.yaml").as_posix()

    with patch(
        "omnibenchmark.cli.run.remote_storage_snakemake_args"
    ) as mock_storage_args:
        with patch(
            "omnibenchmark.cli.run.get_storage_from_benchmark"
        ) as mock_get_storage:
            mock_storage_args.return_value = {"remote-storage": "s3://bucket"}
            mock_get_storage.return_value = MagicMock()

            runner = CliRunner()
            result = runner.invoke(
                run,
                [benchmark_path, "--use-remote-storage", "--yes"],
            )

            # Verify remote storage functions were called
            mock_storage_args.assert_called_once()
            mock_get_storage.assert_called_once()
            mock_workflow_run_workflow.assert_called_once()

            # Verify storage options were passed to workflow
            args, kwargs = mock_workflow_run_workflow.call_args
            assert "remote-storage" in kwargs

            assert result.exit_code == 0


@pytest.mark.short
def test_run_workflow_failure(mock_benchmark_execution, mock_workflow_run_workflow):
    """
    Test that run handles workflow execution failure.
    """
    benchmark_path = Path(data / "mock_benchmark.yaml").as_posix()
    mock_workflow_run_workflow.return_value = False

    runner = CliRunner()
    result = runner.invoke(
        run,
        [benchmark_path, "--yes"],
    )

    # Verify the command failed
    assert result.exit_code == 1


@pytest.mark.short
def test_run_module_with_invalid_timeout():
    """
    Test that run_module handles invalid timeout format correctly.
    """
    benchmark_path = Path(data / "mock_benchmark.yaml").as_posix()

    runner = CliRunner()
    result = runner.invoke(
        run,
        [
            benchmark_path,
            "--module",
            "test_module",
            "--task-timeout",
            "invalid_timeout",
        ],
    )

    # Verify the command failed
    assert result.exit_code == 1


@pytest.mark.short
def test_run_module_dataset_inference_failure(tmp_path):
    """
    Test that run_module fails when dataset cannot be inferred from input files.
    """
    benchmark_path = Path(data / "mock_benchmark.yaml").as_posix()
    test_input_dir = tmp_path / "test_input"
    test_input_dir.mkdir()

    # Create an unrelated file that won't match any dataset
    (test_input_dir / "unrelated_file.txt").write_text("test")

    with patch("omnibenchmark.cli.run.BenchmarkExecution") as mock_execution:
        mock_benchmark = MagicMock()
        mock_node = MagicMock()
        mock_node.is_entrypoint.return_value = False
        mock_node.get_inputs.return_value = ["input.txt"]

        mock_benchmark.get_nodes_by_module_id.return_value = [mock_node]
        mock_benchmark.get_benchmark_datasets.return_value = ["dataset1", "dataset2"]
        mock_execution.return_value = mock_benchmark

        runner = CliRunner()
        result = runner.invoke(
            run,
            [
                benchmark_path,
                "--module",
                "test_module",
                "--input-dir",
                str(test_input_dir),
            ],
        )

        # Verify the command failed with appropriate error
        assert result.exit_code == 1


@pytest.mark.short
def test_run_module_workflow_failure():
    """
    Test that run_module handles workflow execution failure.
    """
    benchmark_path = Path(data / "mock_benchmark.yaml").as_posix()

    with patch("omnibenchmark.cli.run.BenchmarkExecution") as mock_execution:
        with patch(
            "omnibenchmark.cli.run.SnakemakeEngine.run_node_workflow"
        ) as mock_run_node:
            mock_benchmark = MagicMock()
            mock_node = MagicMock()
            mock_node.is_entrypoint.return_value = True
            mock_node.get_inputs.return_value = []

            mock_benchmark.get_nodes_by_module_id.return_value = [mock_node]
            mock_benchmark.get_benchmark_datasets.return_value = ["dataset1"]
            mock_execution.return_value = mock_benchmark

            # Workflow fails
            mock_run_node.return_value = False

            with patch("os.getcwd", return_value="/tmp"):
                with patch("os.listdir", return_value=["dataset1.txt"]):
                    with patch("os.path.exists", return_value=True):
                        with patch("os.path.isdir", return_value=True):
                            runner = CliRunner()
                            result = runner.invoke(
                                run,
                                [
                                    benchmark_path,
                                    "--module",
                                    "dataset1",
                                ],
                            )

                            # Verify the command failed
                            assert result.exit_code == 1


@pytest.mark.short
def test_abort_if_user_does_not_confirm_declined():
    """
    Test that abort_if_user_does_not_confirm raises Abort when user declines.
    """
    from omnibenchmark.cli.run import abort_if_user_does_not_confirm
    from omnibenchmark.cli.utils.logging import logger
    import click

    with patch("click.confirm", return_value=False):
        with pytest.raises(click.Abort):
            abort_if_user_does_not_confirm("test action", logger)


@pytest.mark.short
def test_run_invalid_timeout():
    """
    Test that run handles invalid timeout format.
    """
    benchmark_path = Path(data / "mock_benchmark.yaml").as_posix()

    with patch("omnibenchmark.cli.run.BenchmarkExecution") as mock_execution:
        mock_execution.return_value = MagicMock()

        runner = CliRunner()
        result = runner.invoke(
            run,
            [
                benchmark_path,
                "--task-timeout",
                "not_a_valid_timeout",
                "--yes",
            ],
        )

        # Verify the command failed due to invalid timeout
        assert result.exit_code == 1


@pytest.mark.short
def test_run_module_invalid_input_directory(tmp_path):
    """
    Test that run_module fails when input directory doesn't exist.
    """
    benchmark_path = Path(data / "mock_benchmark.yaml").as_posix()
    nonexistent_dir = tmp_path / "does_not_exist"

    with patch("omnibenchmark.cli.run.BenchmarkExecution") as mock_execution:
        mock_benchmark = MagicMock()
        mock_node = MagicMock()
        mock_node.is_entrypoint.return_value = False

        mock_benchmark.get_nodes_by_module_id.return_value = [mock_node]
        mock_execution.return_value = mock_benchmark

        runner = CliRunner()
        result = runner.invoke(
            run,
            [
                benchmark_path,
                "--module",
                "test_module",
                "--input-dir",
                str(nonexistent_dir),
            ],
        )

        # Click validates the path, so exit code is 2 (usage error)
        assert result.exit_code == 2


@pytest.mark.short
def test_run_module_module_not_found():
    """
    Test that run_module fails when module ID is not found in benchmark.
    """
    benchmark_path = Path(data / "mock_benchmark.yaml").as_posix()

    with patch("omnibenchmark.cli.run.BenchmarkExecution") as mock_execution:
        mock_benchmark = MagicMock()
        # Return empty list - module not found
        mock_benchmark.get_nodes_by_module_id.return_value = []
        mock_execution.return_value = mock_benchmark

        runner = CliRunner()
        result = runner.invoke(
            run,
            [
                benchmark_path,
                "--module",
                "nonexistent_module",
            ],
        )

        # Verify the command failed
        assert result.exit_code == 1


@pytest.mark.short
def test_run_module_missing_required_input_files(tmp_path):
    """
    Test that run_module fails when required input files are missing.
    """
    benchmark_path = Path(data / "mock_benchmark.yaml").as_posix()
    test_input_dir = tmp_path / "test_input"
    test_input_dir.mkdir()

    # Create a dataset file but not the required input
    (test_input_dir / "dataset1.txt").write_text("test")

    with patch("omnibenchmark.cli.run.BenchmarkExecution") as mock_execution:
        mock_benchmark = MagicMock()
        mock_node = MagicMock()
        mock_node.is_entrypoint.return_value = False
        # Require an input file that doesn't exist
        mock_node.get_inputs.return_value = ["{dataset}.required_input.txt"]

        mock_benchmark.get_nodes_by_module_id.return_value = [mock_node]
        mock_benchmark.get_benchmark_datasets.return_value = ["dataset1"]
        mock_execution.return_value = mock_benchmark

        runner = CliRunner()
        result = runner.invoke(
            run,
            [
                benchmark_path,
                "--module",
                "test_module",
                "--input-dir",
                str(test_input_dir),
            ],
        )

        # Verify the command failed due to missing files
        assert result.exit_code == 1


@pytest.mark.short
def test_run_module_entrypoint_without_options():
    """
    Test that run_module works for entrypoint modules without requiring input_dir.
    """
    benchmark_path = Path(data / "mock_benchmark.yaml").as_posix()

    with patch("omnibenchmark.cli.run.BenchmarkExecution") as mock_execution:
        with patch(
            "omnibenchmark.cli.run.SnakemakeEngine.run_node_workflow"
        ) as mock_run_node:
            mock_benchmark = MagicMock()
            mock_node = MagicMock()
            mock_node.is_entrypoint.return_value = True
            mock_node.get_inputs.return_value = []

            mock_benchmark.get_nodes_by_module_id.return_value = [mock_node]
            mock_benchmark.get_benchmark_datasets.return_value = ["dataset1"]
            mock_execution.return_value = mock_benchmark

            mock_run_node.return_value = True

            with patch("os.getcwd", return_value="/tmp"):
                with patch("os.listdir", return_value=["dataset1.txt"]):
                    with patch("os.path.exists", return_value=True):
                        with patch("os.path.isdir", return_value=True):
                            runner = CliRunner()
                            result = runner.invoke(
                                run,
                                [
                                    benchmark_path,
                                    "--module",
                                    "dataset1",
                                ],
                            )

                            # Verify the command succeeded
                            assert result.exit_code == 0
