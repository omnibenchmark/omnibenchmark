import os
from pathlib import Path

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
    expected2 = "Benchmark run has finished successfully."

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
    expected = """Error: Invalid value for 'BENCHMARK_PATH'"""
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
    expected = "Benchmark run has finished successfully"

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
    expected = "Benchmark run has finished successfully"

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
    expected_output = "Running benchmark"
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


def test_local_update_true(tmp_path):
    expected1 = "Running benchmark"
    expected2 = "Benchmark run has finished successfully."
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                str(data / "mock_benchmark.yaml"),
                "--update",
            ],
            input="y",
            cwd=tmp_path,
        )

        assert result.returncode == 0
        assert_in_output(result.stdout, expected1)
        assert_in_output(result.stdout, expected2)


def test_local_update_false():
    # When user declines update, should abort with code 1
    expected = "Are you sure you want to re-run the entire workflow"
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                str(data / "mock_benchmark.yaml"),
                "--update",
            ],
            input="n",
        )

        assert result.returncode == 1
        assert_in_output(result.stdout, expected)


def test_local_dry_update():
    # Dry run with update should complete successfully
    expected = "Running benchmark"
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                str(data / "mock_benchmark.yaml"),
                "--update",
                "--dry",
            ]
        )

        assert result.returncode == 0
        assert_in_output(result.stdout, expected)


def test_benchmark_does_fail_if_one_module_fails(bundled_repos, tmp_path):  # noqa: F811
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                (data / "benchmark_failing_module.yaml").as_posix(),
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
                (data / "benchmark_failing_module.yaml").as_posix(),
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


def test_run_benchmark_with_invalid_timeout():
    """Test that invalid timeout format is rejected."""
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                str(data / "mock_benchmark.yaml"),
                "--task-timeout",
                "invalid_format",
                "--dry",
            ]
        )

        assert result.returncode == 1
        assert (
            "Invalid timeout value" in result.stderr
            or "Invalid timeout value" in result.stdout
        )


def test_run_benchmark_out_dir_with_remote_storage():
    """Test that --out-dir fails when used with --use-remote-storage."""
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                str(data / "mock_benchmark.yaml"),
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
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                str(data / "mock_benchmark.yaml"),
                "--task-timeout",
                "5m",
                "--dry",
            ]
        )

        # Should not fail during timeout parsing
        assert "Invalid timeout value" not in result.stderr
        assert "Invalid timeout value" not in result.stdout
