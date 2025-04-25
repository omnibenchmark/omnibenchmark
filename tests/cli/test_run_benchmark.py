from tests.cli.cli_setup import OmniCLISetup

from .asserts import assert_startswith, assert_in_output
from .fixtures import minio_storage, _minio_container  # noqa: F401
from .path import data


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
                "--local",
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
                "--local",
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
                "--local",
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
                "--local",
            ],
        )

        assert result.returncode == 0
        assert_startswith(result.stdout, expected)


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
                "--local",
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
                "--local",
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
                "--local",
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
                "--local",
                "--update",
                "--dry",
            ]
        )

        assert result.returncode == 0
        assert_startswith(result.stdout, expected)
