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
def test_run_with_double_dash_separator(mock_benchmark_execution):
    """
    Test that arguments after -- are correctly passed to snakemake.
    This verifies the standard execution path (not slurm executor).
    """
    benchmark_path = Path(data / "mock_benchmark.yaml").as_posix()

    with patch("omnibenchmark.cli.run._populate_git_cache"):
        with patch("omnibenchmark.cli.run._generate_explicit_snakefile") as mock_gen:
            mock_gen.return_value = []
            with patch("omnibenchmark.cli.run.write_run_manifest"):
                with patch(
                    "omnibenchmark.cli.run._run_snakemake"
                ) as mock_run_snakemake:
                    mock_run_snakemake.return_value = None

                    mock_benchmark = MagicMock()
                    mock_benchmark.get_benchmark_software_backend.return_value = (
                        MagicMock(value="conda")
                    )
                    mock_benchmark_execution.return_value = mock_benchmark

                    runner = CliRunner()
                    result = runner.invoke(
                        run,
                        [
                            benchmark_path,
                            "--cores",
                            "4",
                            "-k",
                            "--",
                            "--forceall",
                            "--rerun-triggers",
                            "mtime",
                            "--verbose",
                        ],
                    )

                    mock_run_snakemake.assert_called_once()
                    args, kwargs = mock_run_snakemake.call_args

                    assert "extra_snakemake_args" in kwargs
                    extra_args = kwargs["extra_snakemake_args"]
                    assert "--forceall" in extra_args
                    assert "--rerun-triggers" in extra_args
                    assert "mtime" in extra_args
                    assert "--verbose" in extra_args

                    assert kwargs["cores"] == 4
                    assert kwargs["continue_on_error"] is True

                    assert result.exit_code == 0
