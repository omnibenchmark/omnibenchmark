"""Integration tests for git cache with CLI commands."""

import shutil
import tempfile
from pathlib import Path

import pytest

from omnibenchmark.git import describe_cache


@pytest.fixture
def temp_cache_dir(monkeypatch):
    """Create a temporary cache directory for testing."""
    cache_dir = Path(tempfile.mkdtemp()) / "git_cache"
    cache_dir.mkdir(parents=True)

    # Mock the cache directory in all relevant modules
    monkeypatch.setattr("omnibenchmark.git.cache.get_git_cache_dir", lambda: cache_dir)
    monkeypatch.setattr("omnibenchmark.config.get_git_cache_dir", lambda: cache_dir)
    monkeypatch.setattr("omnibenchmark.cli.run.get_git_cache_dir", lambda: cache_dir)

    yield cache_dir

    # Cleanup
    if cache_dir.parent.exists():
        shutil.rmtree(cache_dir.parent)


def test_dry_run_populates_cache(temp_cache_dir):
    """Test that dry-run populates the git cache."""
    from omnibenchmark.benchmark import BenchmarkExecution
    from omnibenchmark.cli.run import _populate_git_cache

    benchmark_path = Path(__file__).parent.parent / "data" / "mock_benchmark.yaml"
    out_dir = Path("out")

    b = BenchmarkExecution(benchmark_path, out_dir)

    # Initially cache should be empty
    initial_status = describe_cache(temp_cache_dir)
    assert len(initial_status) == 0, "Cache should be empty initially"

    # Populate cache
    _populate_git_cache(b)

    # Verify cache was populated
    final_status = describe_cache(temp_cache_dir)
    assert len(final_status) == 2, "Cache should contain 2 repositories"

    # Verify the repositories are correct
    repo_urls = {item["url"] for item in final_status}
    assert "https://github.com/omnibenchmark-example/data.git" in repo_urls
    assert "https://github.com/omnibenchmark-example/process.git" in repo_urls


def test_cache_directory_structure(temp_cache_dir):
    """Test that cache follows Go-style directory structure."""
    from omnibenchmark.benchmark import BenchmarkExecution
    from omnibenchmark.cli.run import _populate_git_cache

    benchmark_path = Path(__file__).parent.parent / "data" / "mock_benchmark.yaml"
    out_dir = Path("out")

    b = BenchmarkExecution(benchmark_path, out_dir)
    _populate_git_cache(b)

    # Check directory structure
    expected_dirs = [
        temp_cache_dir / "github.com" / "omnibenchmark-example" / "data",
        temp_cache_dir / "github.com" / "omnibenchmark-example" / "process",
    ]

    for expected_dir in expected_dirs:
        assert expected_dir.exists(), f"Expected directory {expected_dir} to exist"
        assert expected_dir.is_dir(), f"Expected {expected_dir} to be a directory"

        # Verify it's a git repository
        git_dir = expected_dir / ".git"
        assert git_dir.exists(), f"Expected .git directory in {expected_dir}"
