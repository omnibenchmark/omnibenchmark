import re
import io
import sys
import logging
import tempfile
from pathlib import Path

import pytest

from tests.cli.cli_setup import OmniCLISetup
from tests.io.MinIOStorage_setup import MinIOSetup, TmpMinIOStorage

benchmark_data = Path("..") / "data"
benchmark_data_path = Path(__file__).parent / benchmark_data

minio_testcontainer = MinIOSetup(sys.platform == "linux")
if sys.platform == "linux":
    tempdir = Path(tempfile.gettempdir()) / "ob_test_benchmark004"


@pytest.fixture
def capture_logs():
    """Fixture to capture log output during tests."""
    log_stream = io.StringIO()
    handler = logging.StreamHandler(log_stream)
    logger = logging.getLogger("omnibenchmark")
    logger.addHandler(handler)
    logger.setLevel(logging.DEBUG)

    yield log_stream  # Yield the stream to the test function

    # Cleanup after the test
    logger.removeHandler(handler)
    log_stream.close()


def test_remote(capture_logs):
    expected_output = """
    Benchmark YAML file integrity check passed.
    Running benchmark...
    """
    if not sys.platform == "linux":
        pytest.skip(
            "for GHA, only works on linux (https://docs.github.com/en/actions/using-containerized-services/about-service-containers#about-service-containers)",
            allow_module_level=True,
        )

    with TmpMinIOStorage(minio_testcontainer) as tmp:
        tmp.setup(in_dir=benchmark_data_path, out_dir=tempdir)
        with OmniCLISetup() as omni:
            result = omni.call(
                [
                    "run",
                    "benchmark",
                    "--benchmark",
                    str(tmp.benchmark_file),
                ]
            )

            log_output = capture_logs.getvalue()

            assert result.exit_code == 0
            assert clean(log_output).startswith(clean(expected_output))


def test_benchmark_not_found():
    expected_output = """
    Usage: cli run benchmark [OPTIONS]\nTry 'cli run benchmark --help' for help.\n\nError: Invalid value for '-b' / '--benchmark':
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
        assert result.exit_code == 2
        assert clean(result.output).startswith(clean(expected_output))


def test_benchmark_format_incorrect(capture_logs):
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
        log_output = capture_logs.getvalue()

        assert result.exit_code == 1
        assert clean(log_output) == clean(expected_output)


def test_benchmark_software_does_not_exist(capture_logs):
    expected_output = """
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

        log_output = capture_logs.getvalue()

        assert result.exit_code == 1
        assert clean(log_output) == clean(expected_output)


def test_local(capture_logs):
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

        log_output = capture_logs.getvalue()

        assert result.exit_code == 0
        assert clean(log_output).startswith(clean(expected_output))


def test_local_dry(capture_logs):
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

        log_output = capture_logs.getvalue()

        assert result.exit_code == 0
        assert clean(log_output).startswith(clean(expected_output))


def test_local_update_true(capture_logs):
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

        log_output = capture_logs.getvalue()

        assert result.exit_code == 0
        assert clean(log_output).startswith(clean(expected_output))


def test_local_update_false(capture_logs):
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

        log_output = capture_logs.getvalue()

        assert result.exit_code == 1
        assert clean(log_output).startswith(clean(expected_output))


def test_local_dry_update(capture_logs):
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

        log_output = capture_logs.getvalue()

        assert result.exit_code == 0
        assert clean(log_output).startswith(clean(expected_output))


def clean(output: str) -> str:
    output = output.strip()
    output = output.replace("    ", "")

    # Replace different newline characters with a single '\n'
    normalized_output = re.sub(r"\r\n|\r", "\n", output)

    # Replace multiple spaces and tabs with a single space
    normalized_output = re.sub(r"[ \t]+", " ", normalized_output)

    return normalized_output


def cleanup_buckets_on_exit():
    """Cleanup a testing directory once we are finished."""
    TmpMinIOStorage(minio_testcontainer).cleanup_buckets()


@pytest.fixture(scope="session", autouse=True)
def cleanup(request):
    """Cleanup a testing directory once we are finished."""
    request.addfinalizer(cleanup_buckets_on_exit)
