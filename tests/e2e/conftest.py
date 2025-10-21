import pytest
from tests.fixtures import bundled_repos  # noqa: F401 - pytest fixture


def pytest_addoption(parser):
    # Check if --keep-files option already exists to avoid conflicts
    try:
        parser.addoption(
            "--keep-files",
            action="store_true",
            default=False,
            help="Keep temporary files after test execution for inspection",
        )
    except ValueError:
        # Option already exists, skip adding it
        pass


@pytest.fixture
def e2e_env(request):
    """Fixture that provides e2e test environment parameters."""
    # Check if the request has a param attribute (indirect parametrization)
    if hasattr(request, "param"):
        params = request.param.copy()  # Create a copy to avoid modifying the original
    else:
        # Default values if no parametrization
        params = {"keep_files": False, "current_dir": False}

    # Override with command-line options if provided
    if request.config.getoption("--keep-files"):
        params["keep_files"] = True

    return params


@pytest.fixture
def keep_files(request):
    """Simple fixture to check if files should be kept."""
    return request.config.getoption("--keep-files")
