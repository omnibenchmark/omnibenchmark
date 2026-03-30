import pytest

from tests.remote.S3Storage_setup import RustFSSetup, TmpRustFSStorage

from .path import data


@pytest.fixture
def rustfs_storage(_rustfs_container, tmp_path):
    """Fixture to set up and tear down temporary RustFS storage for each test."""
    # We will use pytest's tmp_path fixture for a unique temporary directory per test
    with TmpRustFSStorage(_rustfs_container) as testcase_storage:
        # Set up a per-test RustFS storage with input and output directories
        testcase_storage.setup(in_dir=data, out_dir=tmp_path)
        yield testcase_storage


# we scope the container fixture to session scope because it is expensive to set up and tear down
@pytest.fixture(scope="session")
def _rustfs_container():
    """Fixture to set up and tear down the RustFS test container for each test."""
    try:
        import docker

        docker.from_env().ping()
    except Exception as e:
        pytest.skip(f"Docker daemon not available: {e}")

    # Initialize a RustFS test container with a lifetime of this test session
    rustfs = RustFSSetup()

    # Yield the container for use in tests
    yield rustfs

    # Here would be a good moment to cleanup after the fixture is used (session scoped).
    # But since this is an ephemeral container, we can save the hassle
    # until we really do need it. Just have this in mind and avoid abusing
    # the test s3 storage for the time being.
