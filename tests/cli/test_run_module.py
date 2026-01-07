from typing import Optional

from .asserts import assert_startswith, clean, assert_in_output
from .cli_setup import OmniCLISetup
from .fixtures import minio_storage, _minio_container  # noqa: F401
from .path import get_benchmark_data_path


def do_first_run(clisetup, file: str, cwd: Optional[str] = None):
    run1 = clisetup.call(["run", file], cwd=cwd)
    assert run1.returncode == 0


def test_run_benchmark_out_dir_with_remote_storage():
    """Test that --out-dir fails when used with --use-remote-storage."""
    path = get_benchmark_data_path()

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                str(path / "mock_benchmark.yaml"),
                "--use-remote-storage",
                "--out-dir",
                "custom_output",
                "--dry",
            ]
        )

        assert result.returncode == 2
        error_msg = "--out-dir can only be used with local storage"
        assert error_msg in result.stderr or error_msg in result.stdout


def test_run_benchmark_with_valid_timeout():
    """Test that valid human-friendly timeout formats are accepted."""
    path = get_benchmark_data_path()

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                str(path / "mock_benchmark.yaml"),
                "--task-timeout",
                "5m",
                "--dry",
            ]
        )

        # Should not fail during timeout parsing
        assert "Invalid timeout value" not in result.stderr
        assert "Invalid timeout value" not in result.stdout


def test_default_entry_module(tmp_path):
    expected = """
        Running module on a local dataset.
        Found 1 workflow nodes for module D1.
        Running module benchmark...
        """

    path = get_benchmark_data_path()

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                str(path / "mock_benchmark.yaml"),
                "--module",
                "D1",
            ],
            cwd=str(tmp_path),
        )

        assert result.returncode == 0
        assert_startswith(result.stdout, expected)


def test_default_nonentry_module_fails(tmp_path):
    """Test running a non-entry module without specifying required options"""
    expected = """
    Running module on a local dataset.
    Found 2 workflow nodes for module P1.
    Error: --input-dir is required for non-entrypoint modules when 'out' folder doesn't exist in current directory.
    """

    path = get_benchmark_data_path()

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                str(path / "mock_benchmark.yaml"),
                "--module",
                "P1",
            ],
            cwd=str(tmp_path),
        )

        assert result.returncode == 1
        assert clean(result.stdout) == clean(expected)


def test_benchmark_not_found():
    path = get_benchmark_data_path()

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                str(path / "does_not_exist.yaml"),
                "--module",
                "D1",
            ]
        )

        # XXX this is too brittle. Should capture exception raised by the module instead.
        expected1 = "Usage: "
        expected2 = "run --help' for help"
        expected3 = "Error: Invalid value for 'BENCHMARK_PATH': Path"
        assert result.returncode == 2
        assert_in_output(result.stderr, expected1)
        assert_in_output(result.stderr, expected2)
        assert_in_output(result.stderr, expected3)


def test_benchmark_format_incorrect():
    path = get_benchmark_data_path()

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                str(path / "benchmark_format_incorrect.yaml"),
                "--module",
                "D1",
            ]
        )
        assert result.returncode == 1


def test_module_not_found():
    expected_output = """
    Running module on a local dataset.
    Error: Could not find module with id `not-existing` in benchmark definition
    """

    path = get_benchmark_data_path()

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                str(path / "mock_benchmark.yaml"),
                "--module",
                "not-existing",
                "--input-dir",
                str(path),
            ]
        )

        assert result.returncode == 1
        assert clean(result.stdout) == clean(expected_output)


def test_behaviour_input(tmp_path):
    expected_output = """
    Running module on a local dataset.
    Found 1 workflow nodes for module D1.
    Running module benchmark...
    """

    path = get_benchmark_data_path()

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                str(path / "mock_benchmark.yaml"),
                "--module",
                "D1",
            ],
            cwd=str(tmp_path),
        )

        assert result.returncode == 0
        assert_startswith(result.stdout, expected_output)


def test_behaviour_input_dry(tmp_path):
    expected_output = """
    Running module on a local dataset.
    Found 1 workflow nodes for module D1.
    Running module benchmark...
    """

    path = get_benchmark_data_path()

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                str(path / "mock_benchmark.yaml"),
                "--module",
                "D1",
                "--dry",
            ],
            cwd=str(tmp_path),
        )

        assert result.returncode == 0
        assert_startswith(result.stdout, expected_output)


def test_behaviour_input_update_true(tmp_path):
    expected_output = """
    Running module on a local dataset.
    Found 1 workflow nodes for module D1.
    Running module benchmark...
    """

    path = get_benchmark_data_path()

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                str(path / "mock_benchmark.yaml"),
                "--module",
                "D1",
                "--update",
            ],
            input="y",
            cwd=str(tmp_path),
        )

        assert result.returncode == 0
        assert_startswith(result.stdout, expected_output)


def test_behaviour_input_update_dry(tmp_path):
    expected_output = """
    Running module on a local dataset.
    Found 1 workflow nodes for module D1.
    Running module benchmark...
    """

    path = get_benchmark_data_path()

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                str(path / "mock_benchmark.yaml"),
                "--module",
                "D1",
                "--update",
                "--dry",
            ],
            input="y",
            cwd=str(tmp_path),
        )

        assert result.returncode == 0
        assert_startswith(result.stdout, expected_output)


def test_behaviour_input_missing_input_dir():
    path = get_benchmark_data_path()

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                str(path / "Benchmark_001.yaml"),
                "--module",
                "M1",
                "--input-dir",
                str(path / "D1" / "default" / "methods"),
            ],
            input="y",
        )

        expected1 = "Usage: "
        expected2 = "run --help' for help"
        expected3 = "Invalid value for '-i' / '--input-dir': Path"

        assert result.returncode == 2
        assert_in_output(result.stderr, expected1)
        assert_in_output(result.stderr, expected2)
        assert_in_output(result.stderr, expected3)


def test_behaviour_input_missing_input_files():
    expected_output = """
    Running module on a local dataset.
    Found 2 workflow nodes for module M1.
    Running module benchmark...
    Error: The following required input files are missing from the input directory: ['D1.meta.json']
    """

    path = get_benchmark_data_path()

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                str(path / "Benchmark_001.yaml"),
                "--module",
                "M1",
                "--input-dir",
                str(path / "D1" / "missing_files"),
            ],
            input="y",
        )

        assert result.returncode == 1
        assert_startswith(result.stdout, expected_output)


def test_behaviour_input_nested_module_dry():
    expected_output = """
    Running module on a local dataset.
    Found 2 workflow nodes for module M1.
    Running module benchmark...
    """

    path = get_benchmark_data_path()

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                str(path / "Benchmark_001.yaml"),
                "--module",
                "M1",
                "--input-dir",
                str(path / "D1" / "default"),
                "--dry",
            ],
            input="y",
        )

        assert result.returncode == 0
        assert_startswith(result.stdout, expected_output)
