from typing import Optional

from tests.cli.test_run_benchmark import assert_in_output

from .asserts import assert_startswith, clean
from .cli_setup import OmniCLISetup
from .fixtures import minio_storage, _minio_container  # noqa: F401
from .path import get_benchmark_data_path


def do_first_run(clisetup, file: str, cwd: Optional[str] = None):
    run1 = clisetup.call(
        [
            "run",
            "benchmark",
            "--benchmark",
            file,
        ],
        cwd=cwd,
    )
    assert run1.returncode == 0


def test_default_entry_module():
    expected = """
        Running module on a local dataset.
        Benchmark YAML file integrity check passed.
        Found 1 workflow nodes for module D1.
        Running module benchmark...
        """

    path = get_benchmark_data_path()

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "module",
                "--benchmark",
                str(path / "mock_benchmark.yaml"),
                "--module",
                "D1",
            ]
        )

        assert result.returncode == 0
        assert_startswith(result.stdout, expected)


def test_default_nonentry_module_fails():
    """Test running a non-entry module without specifying required options"""
    expected = """
    Running module on a local dataset.
    Benchmark YAML file integrity check passed.
    Found 2 workflow nodes for module P1.
    Error: At least one option must be specified. Use '--input_dir', '--example', or '--all'.
    """

    path = get_benchmark_data_path()

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "module",
                "--benchmark",
                str(path / "mock_benchmark.yaml"),
                "--module",
                "P1",
            ]
        )

        assert result.returncode == 1
        assert clean(result.stdout) == clean(expected)


def test_benchmark_not_found():
    path = get_benchmark_data_path()

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "module",
                "--benchmark",
                str(path / "does_not_exist.yaml"),
                "--module",
                "D1",
            ]
        )

        # XXX this is too brittle. Should capture exception raised by the module instead.
        expected1 = "Usage: "
        expected2 = "run module --help' for help"
        expected3 = "Error: Invalid value for '-b' / '--benchmark': Path"
        assert result.returncode == 2
        assert_in_output(result.stderr, expected1)
        assert_in_output(result.stderr, expected2)
        assert_in_output(result.stderr, expected3)


def test_benchmark_format_incorrect():
    expected_output = """
    Running module on a local dataset.
    Error: Failed to parse YAML as a valid OmniBenchmark: software_backend must be supplied.
    """

    path = get_benchmark_data_path()

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "module",
                "--benchmark",
                str(path / "benchmark_format_incorrect.yaml"),
                "--module",
                "D1",
            ]
        )
        assert result.returncode == 1
        assert clean(result.stdout) == clean(expected_output)


def test_module_not_found():
    expected_output = """
    Running module on a local dataset.
    Benchmark YAML file integrity check passed.
    Error: Could not find module with id `not-existing` in benchmark definition
    """

    path = get_benchmark_data_path()

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "module",
                "--benchmark",
                str(path / "mock_benchmark.yaml"),
                "--module",
                "not-existing",
                "--input_dir",
                str(path),
            ]
        )

        assert result.returncode == 1
        assert clean(result.stdout) == clean(expected_output)


def test_behaviour_input():
    expected_output = """
    Running module on a local dataset.
    Benchmark YAML file integrity check passed.
    Found 1 workflow nodes for module D1.
    Running module benchmark...
    """

    path = get_benchmark_data_path()

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "module",
                "--benchmark",
                str(path / "mock_benchmark.yaml"),
                "--module",
                "D1",
            ]
        )

        assert result.returncode == 0
        assert_startswith(result.stdout, expected_output)


def test_behaviour_input_dry():
    expected_output = """
    Running module on a local dataset.
    Benchmark YAML file integrity check passed.
    Found 1 workflow nodes for module D1.
    Running module benchmark...
    """

    path = get_benchmark_data_path()

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "module",
                "--benchmark",
                str(path / "mock_benchmark.yaml"),
                "--module",
                "D1",
                "--dry",
            ]
        )

        assert result.returncode == 0
        assert_startswith(result.stdout, expected_output)


def test_behaviour_input_update_true():
    expected_output = """
    Running module on a local dataset.
    Benchmark YAML file integrity check passed.
    Found 1 workflow nodes for module D1.
    Running module benchmark...
    """

    path = get_benchmark_data_path()

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "module",
                "--benchmark",
                str(path / "mock_benchmark.yaml"),
                "--module",
                "D1",
                "--update",
            ],
            input="y",
        )

        assert result.returncode == 0
        assert_startswith(result.stdout, expected_output)


def test_behaviour_input_update_dry():
    expected_output = """
    Running module on a local dataset.
    Benchmark YAML file integrity check passed.
    Found 1 workflow nodes for module D1.
    Running module benchmark...
    """

    path = get_benchmark_data_path()

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "module",
                "--benchmark",
                str(path / "mock_benchmark.yaml"),
                "--module",
                "D1",
                "--update",
                "--dry",
            ],
            input="y",
        )

        assert result.returncode == 0
        assert_startswith(result.stdout, expected_output)


def test_behaviour_input_missing_input_dir():
    path = get_benchmark_data_path()

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "module",
                "--benchmark",
                str(path / "Benchmark_001.yaml"),
                "--module",
                "M1",
                "--input_dir",
                str(path / "D1" / "default" / "methods"),
            ],
            input="y",
        )

        expected1 = "Usage: "
        expected2 = "run module --help' for help"
        expected3 = "Invalid value for '-i' / '--input_dir': Path"

        assert result.returncode == 2
        assert_in_output(result.stderr, expected1)
        assert_in_output(result.stderr, expected2)
        assert_in_output(result.stderr, expected3)


def test_behaviour_input_missing_input_files():
    expected_output = """
    Running module on a local dataset.
    Benchmark YAML file integrity check passed.
    Found 2 workflow nodes for module M1.
    Running module benchmark...
    Error: The following required input files are missing from the input directory: ['D1.meta.json']
    """

    path = get_benchmark_data_path()

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "module",
                "--benchmark",
                str(path / "Benchmark_001.yaml"),
                "--module",
                "M1",
                "--input_dir",
                str(path / "D1" / "missing_files"),
            ],
            input="y",
        )

        assert result.returncode == 1
        assert_startswith(result.stdout, expected_output)


def test_behaviour_input_nested_module_dry():
    expected_output = """
    Running module on a local dataset.
    Benchmark YAML file integrity check passed.
    Found 2 workflow nodes for module M1.
    Running module benchmark...
    """

    path = get_benchmark_data_path()

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "module",
                "--benchmark",
                str(path / "Benchmark_001.yaml"),
                "--module",
                "M1",
                "--input_dir",
                str(path / "D1" / "default"),
                "--dry",
            ],
            input="y",
        )

        assert result.returncode == 0
        assert_startswith(result.stdout, expected_output)
