import os
import shutil
from pathlib import Path

from tests.cli.cli_setup import OmniCLISetup

from .asserts import assert_startswith, assert_in_output
from .path import data

# TODO: deprecate fixtures in this module
from ..fixtures import minio_storage, _minio_container, bundled_repos  # noqa: F401


# TODO: mark as integration
def test_remote(minio_storage):  # noqa: F811
    # TODO(ben): the technique of expecting YAML validation in the output is a bit brittle, we could
    # check e.g. that output has been produced.
    # But we should be changing the testing strategy in a gradual way
    expected1 = "Benchmark YAML file integrity check passed."
    expected2 = "Benchmark run has finished successfully."

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "benchmark",
                "--benchmark",
                str(minio_storage.benchmark_file),
            ]
        )

        print(result.stdout)
        assert result.returncode == 0
        assert_in_output(result.stdout, expected1)
        assert_in_output(result.stdout, expected2)


def test_benchmark_not_found():
    expected = """Error: Invalid value for '-b' / '--benchmark'"""
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "benchmark",
                "--benchmark",
                str(data / "does_not_exist.yaml"),
                "--local-storage",
            ]
        )
        assert result.returncode == 2
        assert_in_output(result.stderr, expected)


def test_benchmark_format_incorrect():
    # TODO(ben): this should better be a config parsing test
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "benchmark",
                "--benchmark",
                str(data / "benchmark_format_incorrect.yaml"),
                "--local-storage",
            ]
        )
        assert result.returncode == 1


def test_benchmark_software_does_not_exist():
    expected = """
    Error: An unexpected error occurred: Software environment with id 'python' does not have a valid backend definition for: 'conda'.
    """
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "benchmark",
                "--benchmark",
                str(data / "benchmark_software_does_not_exist.yaml"),
                "--local-storage",
            ]
        )

        assert result.returncode == 1
        assert_in_output(result.stdout, expected)


def test_local():
    expected = """
    Benchmark YAML file integrity check passed.
    Running benchmark..."""

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "benchmark",
                "--benchmark",
                str(data / "mock_benchmark.yaml"),
                "--local-storage",
            ],
        )

        assert result.returncode == 0
        assert_startswith(result.stdout, expected)


def test_custom_out_dir():
    expected = """
    Benchmark YAML file integrity check passed.
    Running benchmark..."""

    custom_out_dir = "out_2313_custom"

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "benchmark",
                "--benchmark",
                str(data / "mock_benchmark.yaml"),
                "--local-storage",
                "--out-dir",
                custom_out_dir,
            ],
        )

        assert result.returncode == 0
        assert_startswith(result.stdout, expected)

        assert os.path.exists(custom_out_dir)
        shutil.rmtree(custom_out_dir)


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
                str(data / "mock_benchmark.yaml"),
                "--local-storage",
                "--dry",
            ]
        )

        assert result.returncode == 0
        assert_startswith(result.stdout, expected_output)


def test_local_update_true():
    expected1 = "Running benchmark"
    expected2 = "Benchmark run has finished successfully."
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "benchmark",
                "--benchmark",
                str(data / "mock_benchmark.yaml"),
                "--local-storage",
                "--update",
            ],
            input="y",
        )

        assert result.returncode == 0
        assert_in_output(result.stdout, expected1)
        assert_in_output(result.stdout, expected2)


def test_local_update_false():
    expected = """
    Benchmark YAML file integrity check passed.
    """
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "benchmark",
                "--benchmark",
                str(data / "mock_benchmark.yaml"),
                "--local-storage",
                "--update",
            ],
            input="n",
        )

        assert result.returncode == 1
        assert_startswith(result.stdout, expected)


def test_local_dry_update():
    expected = """
    Benchmark YAML file integrity check passed.
    Running benchmark...
    """
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "benchmark",
                "--benchmark",
                str(data / "mock_benchmark.yaml"),
                "--local-storage",
                "--update",
                "--dry",
            ]
        )

        assert result.returncode == 0
        assert_startswith(result.stdout, expected)


def test_benchmark_does_fail_if_one_module_fails(bundled_repos, tmp_path):  # noqa: F811
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "benchmark",
                "--benchmark",
                (data / "benchmark_failing_module.yaml").as_posix(),
                "--local-storage",
            ],
            input="y",
            cwd=tmp_path,
        )

        failed_msg = "Benchmark run has failed"

        assert failed_msg in result.stdout
        assert result.returncode == 1


def test_benchmark_ok_if_one_module_fails_with_continue(tmp_path, bundled_repos):  # noqa: F811
    """
    This test checks that the benchmark does not fail if one module fails with continue-on-error.
    """
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "benchmark",
                "--benchmark",
                (data / "benchmark_failing_module.yaml").as_posix(),
                "--local-storage",
                "--continue-on-error",
            ],
            input="y",
            cwd=tmp_path,
        )

        failed_msg = "Benchmark run has failed"

        assert failed_msg not in result.stdout
        assert result.returncode == 0

        assert os.path.exists(
            Path(tmp_path)
            / "out"
            / "data/D1/output-D1.txt/process/P1/ok-1/analyze/A1/ok-1_output-analyzed.txt"
        )
