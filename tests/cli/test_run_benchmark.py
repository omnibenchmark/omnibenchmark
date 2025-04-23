import re
import sys
from pathlib import Path


import pytest

from tests.cli.cli_setup import OmniCLISetup
from tests.io.MinIOStorage_setup import MinIOSetup, TmpMinIOStorage


def get_benchmark_data_path() -> Path:
    return Path(__file__).resolve().parent.parent / "data"


benchmark_data_path = get_benchmark_data_path()

# a couple of assertion shortcuts


def assert_startswith(fd, expected):
    assert clean(fd).startswith(clean(expected))


def assert_in_output(fd, expected):
    assert clean(fd) in expected


def clean(output: str) -> str:
    output = output.strip()
    output = output.replace("    ", "")

    # Replace different newline characters with a single '\n'
    normalized_output = re.sub(r"\r\n|\r", "\n", output)

    # Replace multiple spaces and tabs with a single space
    normalized_output = re.sub(r"[ \t]+", " ", normalized_output)

    return normalized_output


@pytest.fixture
def minio_container(scope="session"):
    """Fixture to set up and tear down the MinIO test container for each test."""
    if sys.platform != "linux":
        pytest.skip(
            "for GHA, only works on linux (https://docs.github.com/en/actions/using-containerized-services/about-service-containers#about-service-containers)",
            allow_module_level=True,
        )

    # Initialize a MinIO test container with a lifetime of this test session
    minio = MinIOSetup()

    # Yield the container for use in tests
    yield minio

    # Here would be a good moment to cleanup after the fixture is used (session scoped).
    # But since this is an ephemeral container, we can save the hassle
    # until we really do need it. Just have this in mind and avoid abusing
    # the test s3 storage for the time being.


@pytest.fixture
def minio_storage(minio_container, tmp_path):
    """Fixture to set up and tear down temporary MinIO storage for each test."""
    # We will use pytest's tmp_path fixture for a unique temporary directory per test
    with TmpMinIOStorage(minio_container) as testcase_storage:
        # Set up a per-test MinIO storage with input and output directories
        testcase_storage.setup(in_dir=benchmark_data_path, out_dir=tmp_path)
        yield testcase_storage


# TODO: mark as integration
def test_remote(minio_storage):
    # TODO(ben): the technique of expecting YAML validation in the output is a bit brittle, we could
    # check e.g. that output has been produced.
    # But we should be changing the testing strategy in a gradual way
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
                str(minio_storage.benchmark_file),
            ]
        )

        print(result.stdout)
        assert result.returncode == 0
        assert_startswith(result.stdout, expected)


def test_benchmark_not_found():
    expected = """Error: Invalid value for '-b' / '--benchmark'"""
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
        assert result.returncode == 2
        assert_in_output(result.stdout, expected)


def test_benchmark_format_incorrect():
    expected_output = """
    Error: Failed to parse YAML as a valid OmniBenchmark: software_backend must be supplied.
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

        assert result.returncode == 1
        assert clean(result.stdout) == clean(expected_output)


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
                str(benchmark_data_path / "benchmark_software_does_not_exist.yaml"),
                "--local",
            ]
        )

        assert result.returncode == 1
        assert clean(expected) in result.stdout


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
                str(benchmark_data_path / "mock_benchmark.yaml"),
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
                str(benchmark_data_path / "mock_benchmark.yaml"),
                "--local",
                "--dry",
            ]
        )

        assert result.returncode == 0
        assert_startswith(result.stdout, expected_output)


def test_local_update_true():
    expected_output = """
    Benchmark YAML file integrity check passed.
    Running benchmark...
    Benchmark run has finished successfully.
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
        # TODO: fixme - pass input

        assert result.returncode == 0
        assert_startswith(result.stdout, expected_output)


# TODO: fixme, this is not testing anything really
def test_local_update_false():
    expected_output = """
    Benchmark YAML file integrity check passed.
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

        assert result.returncode == 1
        assert_startswith(result.stdout, expected_output)


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
            ]
        )

        assert result.returncode == 0
        assert_startswith(result.stdout, expected_output)
