import re
from pathlib import Path

from tests.cli.cli_setup import OmniCLISetup

benchmark_path = Path(__file__).parent / ".."
benchmark_data_path = benchmark_path / "data"


def test_default():
    expected_output = """
    Error: At least one option must be specified. Use '--input', '--example', or '--all'.
    """
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "module",
                "--benchmark",
                str(benchmark_data_path / "mock_benchmark.yaml"),
                "--module",
                "D1",
            ]
        )
        assert result.exit_code == 1
        assert clean(result.output) == clean(expected_output)


def test_multiple_behaviours_set():
    expected_output = """
    Error: Only one of '--input', '--example', or '--all' should be set. Please choose only one option.
    """
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "module",
                "--benchmark",
                str(benchmark_data_path / "mock_benchmark.yaml"),
                "--module",
                "D1",
                "--all",
                "--example",
            ]
        )
        assert result.exit_code == 1
        assert clean(result.output) == clean(expected_output)


def test_behaviour_example():
    expected_output = """
    Running module on a predefined remote example dataset.
    Error: Remote execution is not supported yet. Workflows can only be run in local mode.
    """
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "module",
                "--benchmark",
                str(benchmark_data_path / "mock_benchmark.yaml"),
                "--module",
                "D1",
                "--example",
            ]
        )
        assert result.exit_code == 1
        assert clean(result.output) == clean(expected_output)


def test_behaviour_all():
    expected_output = """
    Running module on all available remote datasets.
    Error: Remote execution is not supported yet. Workflows can only be run in local mode.
    """
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "module",
                "--benchmark",
                str(benchmark_data_path / "mock_benchmark.yaml"),
                "--module",
                "D1",
                "--all",
            ]
        )
        assert result.exit_code == 1
        assert clean(result.output) == clean(expected_output)


def test_benchmark_not_found():
    expected_output = """
    Running module on a dataset provided in a custom directory.
    Error: Benchmark YAML file not found.
    """
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "module",
                "--benchmark",
                str(benchmark_data_path / "does_not_exist.yaml"),
                "--module",
                "D1",
                "--input",
                str(benchmark_path),
            ]
        )
        assert result.exit_code == 1
        assert clean(result.output) == clean(expected_output)


def test_benchmark_format_incorrect():
    expected_output = """
    Running module on a dataset provided in a custom directory.
    Error: Failed to parse YAML as a valid OmniBenchmark: software_environments must be supplied.
    """
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "module",
                "--benchmark",
                str(benchmark_data_path / "benchmark_format_incorrect.yaml"),
                "--module",
                "D1",
                "--input",
                str(benchmark_path),
            ]
        )
        assert result.exit_code == 1
        assert clean(result.output) == clean(expected_output)


def test_module_not_found():
    expected_output = """
    Running module on a dataset provided in a custom directory.
    Benchmark YAML file integrity check passed.
    Error: Could not find module with id `not-existing` in benchmark definition
    """
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "module",
                "--benchmark",
                str(benchmark_data_path / "mock_benchmark.yaml"),
                "--module",
                "not-existing",
                "--input",
                str(benchmark_path),
            ]
        )
        assert result.exit_code == 1
        assert clean(result.output) == clean(expected_output)


def test_behaviour_input():
    expected_output = """
    Running module on a dataset provided in a custom directory.
    Benchmark YAML file integrity check passed.
    Found 1 workflow nodes for module D1.
    Running module benchmark...
    """
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "module",
                "--benchmark",
                str(benchmark_data_path / "mock_benchmark.yaml"),
                "--module",
                "D1",
                "--input",
                str(benchmark_path),
            ]
        )
        assert result.exit_code == 0
        assert clean(result.output).startswith(clean(expected_output))


def test_behaviour_input_dry():
    expected_output = """
    Running module on a dataset provided in a custom directory.
    Benchmark YAML file integrity check passed.
    Found 1 workflow nodes for module D1.
    Running module benchmark...
    """
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "module",
                "--benchmark",
                str(benchmark_data_path / "mock_benchmark.yaml"),
                "--module",
                "D1",
                "--input",
                str(benchmark_path),
                "--dry",
            ]
        )
        assert result.exit_code == 0
        assert clean(result.output).startswith(clean(expected_output))


def test_behaviour_input_update_true():
    expected_output = """
    Running module on a dataset provided in a custom directory.
    Benchmark YAML file integrity check passed.
    Are you sure you want to re-run the entire workflow? [y/N]: y
    Found 1 workflow nodes for module D1.
    Running module benchmark...
    """
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "module",
                "--benchmark",
                str(benchmark_data_path / "mock_benchmark.yaml"),
                "--module",
                "D1",
                "--input",
                str(benchmark_path),
                "--update",
            ],
            input="y",
        )
        assert result.exit_code == 0
        assert clean(result.output).startswith(clean(expected_output))


def test_behaviour_input_update_false():
    expected_output = """
    Running module on a dataset provided in a custom directory.
    Benchmark YAML file integrity check passed.
    Are you sure you want to re-run the entire workflow? [y/N]: N
    """
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "module",
                "--benchmark",
                str(benchmark_data_path / "mock_benchmark.yaml"),
                "--module",
                "D1",
                "--input",
                str(benchmark_path),
                "--update",
            ],
            input="N",
        )
        assert result.exit_code == 1
        assert clean(result.output).startswith(clean(expected_output))


def test_behaviour_input_update_dry():
    expected_output = """
    Running module on a dataset provided in a custom directory.
    Benchmark YAML file integrity check passed.
    Found 1 workflow nodes for module D1.
    Running module benchmark...
    """
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "module",
                "--benchmark",
                str(benchmark_data_path / "mock_benchmark.yaml"),
                "--module",
                "D1",
                "--input",
                str(benchmark_path),
                "--update",
                "--dry",
            ],
            input="y",
        )
        assert result.exit_code == 0
        assert clean(result.output).startswith(clean(expected_output))


def test_behaviour_input_missing_input_dir():
    expected_output = """
    Running module on a dataset provided in a custom directory.
    Benchmark YAML file integrity check passed.
    Found 2 workflow nodes for module M1.
    Running module benchmark...
    Error: Input directory does not exist or is not a valid directory:
    """
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "module",
                "--benchmark",
                str(benchmark_data_path / "Benchmark_001.yaml"),
                "--module",
                "M1",
                "--input",
                str(benchmark_path / "data" / "D1" / "default" / "methods"),
            ],
            input="y",
        )
        assert result.exit_code == 1
        assert clean(result.output).startswith(clean(expected_output))


def test_behaviour_input_missing_input_files():
    expected_output = """
    Running module on a dataset provided in a custom directory.
    Benchmark YAML file integrity check passed.
    Found 2 workflow nodes for module M1.
    Running module benchmark...
    Error: The following required input files are missing from the input directory: ['D1.meta.json']
    """
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "module",
                "--benchmark",
                str(benchmark_data_path / "Benchmark_001.yaml"),
                "--module",
                "M1",
                "--input",
                str(benchmark_path / "data" / "D1" / "missing_files"),
            ],
            input="y",
        )
        assert result.exit_code == 1
        assert clean(result.output).startswith(clean(expected_output))


def test_behaviour_input_nested_module_dry():
    expected_output = """
    Running module on a dataset provided in a custom directory.
    Benchmark YAML file integrity check passed.
    Found 2 workflow nodes for module M1.
    Running module benchmark...
    """
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "module",
                "--benchmark",
                str(benchmark_data_path / "Benchmark_001.yaml"),
                "--module",
                "M1",
                "--input",
                str(benchmark_path / "data" / "D1" / "default"),
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
