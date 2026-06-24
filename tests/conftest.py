import io
import os
import shutil
import tempfile

# Isolate the global omnibenchmark git cache from the developer's real cache and
# from stale cross-run state.  The cache (normally ~/.cache/omnibenchmark/git) is
# keyed partly by ephemeral tmp paths, so entries left by previous runs point at
# deleted tmp dirs and can make `git` updates fail intermittently — surfacing as
# flaky cli/e2e failures.  Redirecting XDG_CACHE_HOME to a fresh per-session dir
# *before* omnibenchmark is imported (config freezes the cache path at import
# time) fixes it for in-process code and for the CLI subprocesses that inherit
# this environment.  Cleaned up in pytest_unconfigure.
_TEST_XDG_CACHE_HOME = tempfile.mkdtemp(prefix="ob-test-xdg-cache-")
os.environ["XDG_CACHE_HOME"] = _TEST_XDG_CACHE_HOME

import pytest  # noqa: E402
import logging  # noqa: E402

from pathlib import Path  # noqa: E402
from .git_bundle import GitBundleManager  # noqa: E402


def pytest_unconfigure(config):
    """Remove the isolated git cache dir created at import time."""
    shutil.rmtree(_TEST_XDG_CACHE_HOME, ignore_errors=True)


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
