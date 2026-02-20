import io
import os

import pytest
import logging

from pathlib import Path
from .git_bundle import GitBundleManager


@pytest.fixture(autouse=True)
def isolated_tmp_cwd(tmp_path, request):
    """Run every test in an isolated temporary working directory.

    Prevents stray output files (out/, .modules/, etc.) in the project root
    from leaking into tests that read the filesystem relative to cwd.

    Tests that explicitly need the project root (e.g. e2e tests that reference
    fixture paths) can opt out with @pytest.mark.no_tmp_cwd.
    """
    if request.node.get_closest_marker("no_tmp_cwd"):
        yield
        return
    original = Path.cwd()
    os.chdir(tmp_path)
    try:
        yield tmp_path
    finally:
        os.chdir(original)


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


# git fixtures


@pytest.fixture(scope="session")
def bundles_dir() -> Path:
    """Directory containing test bundles."""
    return Path(__file__).parent / "data" / "bundles"


@pytest.fixture(scope="session")
def repos_dir(tmp_path_factory) -> Path:
    """Directory where repositories will be extracted."""
    return Path(__file__).parent / "data" / "repos"


@pytest.fixture(scope="session")
def git_bundle_manager(bundles_dir, repos_dir):
    """Fixture providing Git bundle management."""
    return GitBundleManager(bundles_dir, repos_dir)


@pytest.fixture
def data_repo(git_bundle_manager):
    """Fixture providing an extracted data repository."""
    repo_path = git_bundle_manager.extract_bundle("data")
    if repo_path is None:
        pytest.skip("Data bundle not available")
    return repo_path
