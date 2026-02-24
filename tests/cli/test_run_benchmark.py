import os

import pytest

from tests.cli.cli_setup import OmniCLISetup

from .asserts import assert_in_output
from .path import data

# TODO: deprecate fixtures in this module
from ..fixtures import minio_storage, _minio_container, bundled_repos  # noqa: F401


# TODO: mark as integration
def test_remote(minio_storage):  # noqa: F811
    # TODO(ben): the technique of expecting YAML validation in the output is a bit brittle, we could
    # check e.g. that output has been produced.
    # But we should be changing the testing strategy in a gradual way
    expected1 = "Benchmark YAML file integrity check passed."
    expected2 = "Benchmark run completed successfully."

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                str(minio_storage.benchmark_file),
            ]
        )

        print(result.stdout)
        assert result.returncode == 0
        assert_in_output(result.stdout, expected1)
        assert_in_output(result.stdout, expected2)


def test_benchmark_not_found():
    expected = """Error: Invalid value for 'BENCHMARK'"""
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                str(data / "does_not_exist.yaml"),
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
                str(data / "benchmark_format_incorrect.yaml"),
            ]
        )
        assert result.returncode == 1


def test_benchmark_software_does_not_exist():
    expected = """
    Failed to load benchmark: Software environment with id 'python' does not have a valid backend definition for: 'conda'.
    """
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                str(data / "benchmark_software_does_not_exist.yaml"),
            ]
        )

        assert result.returncode == 1
        assert_in_output(result.stdout, expected)


def test_local(tmp_path):
    # Check that benchmark runs successfully (may have deprecation warnings)
    expected = "Benchmark run completed successfully."

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                str(data / "mock_benchmark.yaml"),
            ],
            cwd=tmp_path,
        )

        assert result.returncode == 0
        assert_in_output(result.stdout, expected)


def test_custom_out_dir(tmp_path):
    # Check that benchmark runs successfully with custom output directory
    expected = "Benchmark run completed successfully."

    custom_out_dir = "out_2313_custom"

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                str(data / "mock_benchmark.yaml"),
                "--out-dir",
                custom_out_dir,
            ],
            cwd=tmp_path,
        )

        assert result.returncode == 0
        assert_in_output(result.stdout, expected)

        assert os.path.exists(tmp_path / custom_out_dir)


def test_local_dry():
    # Dry run should complete successfully (may have warnings)
    expected_output = "Snakefile generated."
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                str(data / "mock_benchmark.yaml"),
                "--dry",
            ]
        )

        assert result.returncode == 0
        assert_in_output(result.stdout, expected_output)


@pytest.mark.skip(
    reason="bundle-based repos not supported by the new explicit-snakefile resolver; "
    "requires migration of test data to real git repos"
)
def test_benchmark_does_fail_if_one_module_fails(bundled_repos, tmp_path):  # noqa: F811
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                (data / "benchmark_failing_module.yaml").as_posix(),
            ],
            cwd=tmp_path,
        )

        failed_msg = "Benchmark run failed"

        assert failed_msg in result.stdout or failed_msg in result.stderr
        assert result.returncode != 0


@pytest.mark.skip(
    reason="bundle-based repos not supported by the new explicit-snakefile resolver; "
    "requires migration of test data to real git repos"
)
def test_benchmark_ok_if_one_module_fails_with_continue(tmp_path, bundled_repos):  # noqa: F811
    """
    This test checks that the benchmark does not fail if one module fails with continue-on-error.
    """
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                (data / "benchmark_failing_module.yaml").as_posix(),
                "--continue-on-error",
            ],
            cwd=tmp_path,
        )

        failed_msg = "Benchmark run failed"

        assert failed_msg not in result.stdout
        assert failed_msg not in result.stderr
        assert result.returncode == 0
