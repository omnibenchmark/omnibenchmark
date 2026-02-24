from typing import Optional

import pytest

from .asserts import assert_in_output
from .cli_setup import OmniCLISetup
from .fixtures import minio_storage, _minio_container  # noqa: F401
from .path import get_benchmark_data_path


def do_first_run(clisetup, file: str, cwd: Optional[str] = None):
    run1 = clisetup.call(["run", file], cwd=cwd)
    assert run1.returncode == 0


@pytest.mark.skip(reason="--use-remote-storage removed from new explicit-snakefile CLI")
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
    """Test running an entry module via --module filter."""
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
        assert "Snakefile generated." in result.stdout or "Generated" in result.stdout


@pytest.mark.skip(
    reason="--input-dir option removed from new explicit-snakefile CLI; "
    "module filter no longer checks for non-entrypoint requirements"
)
def test_default_nonentry_module_fails(tmp_path):
    """Test running a non-entry module without specifying required options"""
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
        expected3 = "Error: Invalid value for 'BENCHMARK': Path"
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
    path = get_benchmark_data_path()

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                str(path / "mock_benchmark.yaml"),
                "--module",
                "not-existing",
            ]
        )

        assert result.returncode == 1
        assert "not-existing" in result.stdout or "not-existing" in result.stderr


def test_behaviour_input(tmp_path):
    """Test running benchmark with --module filter (entry-point module)."""
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


def test_behaviour_input_dry(tmp_path):
    """Test dry run with --module filter."""
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
        assert "Snakefile generated." in result.stdout or "Generated" in result.stdout


@pytest.mark.skip(reason="--update removed from new explicit-snakefile CLI")
def test_behaviour_input_update_true(tmp_path):
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


@pytest.mark.skip(reason="--update removed from new explicit-snakefile CLI")
def test_behaviour_input_update_dry(tmp_path):
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


@pytest.mark.skip(reason="--input-dir removed from new explicit-snakefile CLI")
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

        assert result.returncode == 2


@pytest.mark.skip(reason="--input-dir removed from new explicit-snakefile CLI")
def test_behaviour_input_missing_input_files():
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


@pytest.mark.skip(reason="--input-dir removed from new explicit-snakefile CLI")
def test_behaviour_input_nested_module_dry():
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
