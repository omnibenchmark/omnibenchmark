import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--keep-files",
        action="store_true",
        default=False,
        help="Keep temporary files after test execution",
    )
    parser.addoption(
        "--current-dir",
        action="store_true",
        default=False,
        help="Use current directory instead of temporary one",
    )


@pytest.fixture
def snakemake_env(request):
    """Fixture that provides snakemake environment parameters.
    Can be overridden by command-line options."""
    # Check if the request has a param attribute (indirect parametrization)
    if hasattr(request, "param"):
        params = request.param.copy()  # Create a copy to avoid modifying the original
    else:
        # Default values if no parametrization
        params = {"benchmark_file": None, "keep_files": False, "current_dir": False}

    # Override with command-line options if provided
    if request.config.getoption("--keep-files"):
        params["keep_files"] = True

    if request.config.getoption("--current-dir"):
        params["current_dir"] = True

    return params
