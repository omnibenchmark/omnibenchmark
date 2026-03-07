"""Tests for the new git caching system."""

import pytest
from dulwich import porcelain
from omnibenchmark.git.cache import (
    copy_local_to_work_dir,
    is_local_path,
    parse_repo_url,
    get_or_update_cached_repo,
    checkout_to_work_dir,
    clone_module_v2,
    resolve_local_path,
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


@pytest.fixture
def local_git_repo(tmp_path):
    """Create a minimal local git repo with a commit."""
    repo_dir = tmp_path / "myrepo"
    repo_dir.mkdir()
    porcelain.init(str(repo_dir))
    # Add a file and commit
    (repo_dir / "hello.txt").write_text("hello")
    porcelain.add(str(repo_dir), paths=["hello.txt"])
    commit_sha = porcelain.commit(
        str(repo_dir),
        message=b"initial commit",
        author=b"Test <test@test>",
        committer=b"Test <test@test>",
    )
    return repo_dir, commit_sha.decode("ascii")


class TestIsLocalPath:
    """Test detection of local filesystem paths."""

    def test_dot(self):
        assert is_local_path(".") is True

    def test_dotdot(self):
        assert is_local_path("..") is True

    def test_dot_slash(self):
        assert is_local_path("./subdir") is True

    def test_dotdot_slash(self):
        assert is_local_path("../sibling") is True

    def test_absolute_path(self):
        assert is_local_path("/home/user/repo") is True

    def test_file_url(self):
        assert is_local_path("file:///home/user/repo") is True

    def test_https_url(self):
        assert is_local_path("https://github.com/user/repo") is False

    def test_ssh_url(self):
        assert is_local_path("git@github.com:user/repo") is False

    def test_whitespace_stripped(self):
        assert is_local_path("  .  ") is True


class TestResolveLocalPath:
    """Test resolving local paths to absolute paths."""

    def test_dot_resolves_to_benchmark_dir(self, tmp_path):
        result = resolve_local_path(".", benchmark_dir=tmp_path)
        assert result == tmp_path.resolve()

    def test_relative_subdir(self, tmp_path):
        subdir = tmp_path / "sub"
        subdir.mkdir()
        result = resolve_local_path("./sub", benchmark_dir=tmp_path)
        assert result == subdir.resolve()

    def test_absolute_path_unchanged(self, tmp_path):
        result = resolve_local_path(str(tmp_path), benchmark_dir=None)
        assert result == tmp_path.resolve()

    def test_file_url_stripped(self, tmp_path):
        result = resolve_local_path(f"file://{tmp_path}", benchmark_dir=None)
        assert result == tmp_path.resolve()


class TestCopyLocalToWorkDir:
    """Test copying a local repo to a work directory."""

    def test_copies_files_without_git(self, local_git_repo, tmp_path):
        repo_dir, commit = local_git_repo
        work_dir = tmp_path / "work" / "out"

        result_dir, result_commit = copy_local_to_work_dir(repo_dir, work_dir)

        assert result_dir == work_dir
        assert (work_dir / "hello.txt").exists()
        assert (work_dir / "hello.txt").read_text() == "hello"
        assert not (work_dir / ".git").exists()
        assert result_commit == commit

    def test_nonexistent_path_raises(self, tmp_path):
        with pytest.raises(RuntimeError, match="does not exist"):
            copy_local_to_work_dir(tmp_path / "nope", tmp_path / "work")

    def test_non_git_dir_raises(self, tmp_path):
        plain_dir = tmp_path / "plain"
        plain_dir.mkdir()
        with pytest.raises(RuntimeError, match="not a valid git repository"):
            copy_local_to_work_dir(plain_dir, tmp_path / "work")

    def test_overwrites_existing_work_dir(self, local_git_repo, tmp_path):
        repo_dir, _ = local_git_repo
        work_dir = tmp_path / "work" / "out"
        work_dir.mkdir(parents=True)
        (work_dir / "stale.txt").write_text("stale")

        copy_local_to_work_dir(repo_dir, work_dir)

        assert not (work_dir / "stale.txt").exists()
        assert (work_dir / "hello.txt").exists()
