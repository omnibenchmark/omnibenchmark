import re
from pathlib import Path

from tests.cli.cli_setup import OmniCLISetup

benchmark_data = Path("..") / "data"
benchmark_data_path = Path(__file__).parent / benchmark_data


def test_default():
    expected_output = """
    Error: Remote execution is not supported yet. Workflows can only be run in local mode.
    """
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "benchmark",
                "--benchmark",
                str(benchmark_data_path / "mock_benchmark.yaml"),
            ]
        )
        assert result.exit_code == 1
        assert clean(result.output) == clean(expected_output)


def test_benchmark_not_found():
    expected_output = """
    Error: Benchmark YAML file not found.
    """
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "benchmark",
                "--benchmark",
                str(benchmark_data_path / "does_not_exist.yaml"),
                "--local",
            ]
        )
        assert result.exit_code == 1
        assert clean(result.output) == clean(expected_output)


def test_benchmark_format_incorrect():
    expected_output = """
    Error: Failed to parse YAML as a valid OmniBenchmark: software_environments must be supplied.
    """
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "benchmark",
                "--benchmark",
                str(benchmark_data_path / "benchmark_format_incorrect.yaml"),
                "--local",
            ]
        )
        assert result.exit_code == 1
        assert clean(result.output) == clean(expected_output)


def test_local():
    expected_output = """
    Benchmark YAML file integrity check passed.
    Running benchmark...
    """
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "benchmark",
                "--benchmark",
                str(benchmark_data_path / "mock_benchmark.yaml"),
                "--local",
            ],
        )
        assert result.exit_code == 0
        assert clean(result.output).startswith(clean(expected_output))


def test_local_dry():
    expected_output = """
    Benchmark YAML file integrity check passed.
    Running benchmark...
    """
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "benchmark",
                "--benchmark",
                str(benchmark_data_path / "mock_benchmark.yaml"),
                "--local",
                "--dry",
            ]
        )
        assert result.exit_code == 0
        assert clean(result.output).startswith(clean(expected_output))


def test_local_update_true():
    expected_output = """
    Benchmark YAML file integrity check passed.
    Are you sure you want to re-run the entire workflow? [y/N]: Y
    Running benchmark...
    """
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "benchmark",
                "--benchmark",
                str(benchmark_data_path / "mock_benchmark.yaml"),
                "--local",
                "--update",
            ],
            input="Y",
        )
        assert result.exit_code == 0
        assert clean(result.output).startswith(clean(expected_output))


def test_local_update_false():
    expected_output = """
    Benchmark YAML file integrity check passed.
    Are you sure you want to re-run the entire workflow? [y/N]: n
    """
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "benchmark",
                "--benchmark",
                str(benchmark_data_path / "mock_benchmark.yaml"),
                "--local",
                "--update",
            ],
            input="n",
        )
        assert result.exit_code == 1
        assert clean(result.output).startswith(clean(expected_output))


def test_local_dry_update():
    expected_output = """
    Benchmark YAML file integrity check passed.
    Running benchmark...
    """
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "benchmark",
                "--benchmark",
                str(benchmark_data_path / "mock_benchmark.yaml"),
                "--local",
                "--update",
                "--dry",
            ],
            input="y",
        )
        assert result.exit_code == 0
        assert clean(result.output).startswith(clean(expected_output))


def clean(output: str) -> str:
    output = output.strip()
    output = output.replace("    ", "")

    # Replace different newline characters with a single '\n'
    normalized_output = re.sub(r"\r\n|\r", "\n", output)

    # Replace multiple spaces and tabs with a single space
    normalized_output = re.sub(r"[ \t]+", " ", normalized_output)

    return normalized_output
