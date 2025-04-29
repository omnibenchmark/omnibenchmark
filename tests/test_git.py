# Test git repos support functionality

import subprocess

import pytest


@pytest.fixture
def data_repo(git_bundle_manager):
    return git_bundle_manager.extract_bundle("data_63b7b36")


def test_data_repo_contents(data_repo):
    """Test extractaction of bundle data repository."""
    config = data_repo / "config.cfg"
    assert config.exists(), "config.cfg missing"

    requirements = data_repo / "requirements.txt"
    assert requirements.exists(), "requirements.txt missing"

    entrypoint = data_repo / "entrypoint_data.py"
    assert entrypoint.exists(), "entrypoint_data.py missing"

    # Verify we're at the correct commit
    result = subprocess.run(
        ["git", "rev-parse", "--short=7", "HEAD"],
        cwd=data_repo,
        capture_output=True,
        text=True,
        check=True,
    )
    assert "63b7b36" in result.stdout, "Wrong commit checked out"
