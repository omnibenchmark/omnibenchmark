import sys

import pytest

from pathlib import Path

from tests.remote.MinIOStorage_setup import MinIOSetup, TmpMinIOStorage


def get_benchmark_data_path() -> Path:
    return Path(__file__).resolve().parent / "data"


data = get_benchmark_data_path()


@pytest.fixture
def minio_storage(_minio_container, tmp_path):
    """Fixture to set up and tear down temporary MinIO storage for each test."""
    # We will use pytest's tmp_path fixture for a unique temporary directory per test
    with TmpMinIOStorage(_minio_container) as testcase_storage:
        # Set up a per-test MinIO storage with input and output directories
        testcase_storage.setup(in_dir=data, out_dir=tmp_path)
        yield testcase_storage


# we scope the minio container fixture to session scope because it is expensive to set up and tear down
@pytest.fixture(scope="session")
def _minio_container():
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
def bundled_repos(tmp_path):
    """
    Creates a symlink from tmp_path/bundles to tests/data/bundles.

    This allows tests to access predefined bundle files through a path that appears
    to be in the temporary directory used for each test.
    """
    # Get the source directory (tests/data/bundles)
    source_dir = data / "bundles"

    # Create the target directory in tmp_path
    target_dir = tmp_path / "bundles"

    # Create the symlink
    target_dir.symlink_to(source_dir, target_is_directory=True)

    return target_dir
