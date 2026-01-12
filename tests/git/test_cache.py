"""Tests for the new git caching system."""

import pytest
from omnibenchmark.git.cache import (
    parse_repo_url,
    get_or_update_cached_repo,
    checkout_to_work_dir,
    clone_module_v2,
)


class TestParseRepoUrl:
    """Test URL parsing for Go-style cache paths."""

    def test_https_github(self):
        url = "https://github.com/user/repo.git"
        assert parse_repo_url(url) == "github.com/user/repo"

    def test_https_github_no_git_suffix(self):
        url = "https://github.com/user/repo"
        assert parse_repo_url(url) == "github.com/user/repo"

    def test_ssh_github(self):
        url = "git@github.com:user/repo.git"
        assert parse_repo_url(url) == "github.com/user/repo"

    def test_ssh_github_no_git_suffix(self):
        url = "git@github.com:user/repo"
        assert parse_repo_url(url) == "github.com/user/repo"

    def test_gitlab_nested(self):
        url = "https://gitlab.com/group/subgroup/project.git"
        assert parse_repo_url(url) == "gitlab.com/group/subgroup/project"

    def test_trailing_slash(self):
        url = "https://github.com/user/repo/"
        assert parse_repo_url(url) == "github.com/user/repo"


class TestGetOrUpdateCachedRepo:
    """Test caching full repositories."""

    @pytest.mark.integration
    def test_cache_new_repo(self, tmp_path):
        """Test caching a new repository."""
        repo_url = "https://github.com/omnibenchmark/clustering_example.git"
        cache_dir = tmp_path / "cache"

        repo_path = get_or_update_cached_repo(repo_url, cache_dir)

        # Verify structure: cache/github.com/omnibenchmark/clustering_example
        assert repo_path.exists()
        assert "github.com/omnibenchmark/clustering_example" in str(repo_path)
        assert (repo_path / ".git").exists()

    @pytest.mark.integration
    def test_cache_idempotent(self, tmp_path):
        """Test that caching the same repo twice reuses the cache."""
        repo_url = "https://github.com/omnibenchmark/clustering_example.git"
        cache_dir = tmp_path / "cache"

        # First cache
        path1 = get_or_update_cached_repo(repo_url, cache_dir)

        # Second cache should return same path and fetch updates
        path2 = get_or_update_cached_repo(repo_url, cache_dir)

        assert path1 == path2
        assert path1.exists()


class TestCheckoutToWorkDir:
    """Test checking out to work directories."""

    @pytest.mark.integration
    def test_checkout_by_branch(self, tmp_path):
        """Test checking out a branch."""
        repo_url = "https://github.com/omnibenchmark/clustering_example.git"
        ref = "main"
        cache_dir = tmp_path / "cache"
        work_dir = tmp_path / "work" / "test1"

        result_path, commit = checkout_to_work_dir(repo_url, ref, work_dir, cache_dir)

        assert result_path == work_dir
        assert work_dir.exists()
        assert len(commit) == 40
        assert all(c in "0123456789abcdef" for c in commit)

    @pytest.mark.integration
    def test_checkout_by_commit(self, tmp_path):
        """Test checking out a specific commit."""
        repo_url = "https://github.com/omnibenchmark/clustering_example.git"
        commit = "1faafa2"  # Short hash
        cache_dir = tmp_path / "cache"
        work_dir = tmp_path / "work" / "test2"

        result_path, resolved_commit = checkout_to_work_dir(
            repo_url, commit, work_dir, cache_dir
        )

        assert result_path == work_dir
        assert work_dir.exists()
        assert resolved_commit.lower().startswith(commit.lower())


class TestCloneModuleV2:
    """Test the main entry point for module cloning."""

    @pytest.mark.integration
    def test_clone_by_branch(self, tmp_path):
        """Test cloning by branch name."""
        repo_url = "https://github.com/omnibenchmark/clustering_example.git"
        branch = "main"
        work_dir = tmp_path / "work" / "branch_test"
        cache_dir = tmp_path / "cache"

        module_path, commit = clone_module_v2(repo_url, branch, work_dir, cache_dir)

        assert module_path == work_dir
        assert module_path.exists()
        assert len(commit) == 40
        assert all(c in "0123456789abcdef" for c in commit)

    @pytest.mark.integration
    def test_clone_by_commit(self, tmp_path):
        """Test cloning by commit hash."""
        repo_url = "https://github.com/omnibenchmark/clustering_example.git"
        commit = "1faafa2"
        work_dir = tmp_path / "work" / "commit_test"
        cache_dir = tmp_path / "cache"

        module_path, resolved_commit = clone_module_v2(
            repo_url, commit, work_dir, cache_dir
        )

        assert module_path == work_dir
        assert module_path.exists()
        assert resolved_commit.lower().startswith(commit.lower())

    @pytest.mark.integration
    def test_go_style_cache_structure(self, tmp_path):
        """Verify the Go-style cache structure."""
        repo_url = "https://github.com/omnibenchmark/clustering_example.git"
        commit = "1faafa2"
        work_dir = tmp_path / "work" / "structure_test"
        cache_dir = tmp_path / "cache"

        module_path, _ = clone_module_v2(repo_url, commit, work_dir, cache_dir)

        # Verify cache structure: cache/github.com/omnibenchmark/clustering_example
        cache_repo_path = (
            cache_dir / "github.com" / "omnibenchmark" / "clustering_example"
        )
        assert cache_repo_path.exists()
        assert (cache_repo_path / ".git").exists()

        # Verify work dir is separate
        assert module_path.exists()
        assert module_path != cache_repo_path
