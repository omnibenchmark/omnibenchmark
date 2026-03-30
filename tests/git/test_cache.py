"""Tests for the new git caching system."""

import random
import threading
from unittest.mock import patch

import pytest
from dulwich import porcelain
from dulwich.repo import Repo

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


# ---------------------------------------------------------------------------
# Fixtures for checkout_to_work_dir concurrency tests
# ---------------------------------------------------------------------------


@pytest.fixture
def two_commit_cache(tmp_path):
    """
    A local git repo with two commits (alpha, beta) and a pre-cloned cache copy.

    Returns (cache_path, sha_alpha, sha_beta) where:
      - sha_alpha has content.txt = "version: alpha\\n" and extra.txt
      - sha_beta  has content.txt = "version: beta\\n"  and extra.txt
    """
    # Build origin repo with two sequential commits
    origin = tmp_path / "origin"
    origin.mkdir()
    porcelain.init(str(origin))

    (origin / "content.txt").write_text("version: alpha\n")
    (origin / "extra.txt").write_text("extra file\n")
    porcelain.add(str(origin), paths=["content.txt", "extra.txt"])
    sha_alpha = porcelain.commit(
        str(origin),
        message=b"alpha commit",
        author=b"Test <test@example.com>",
        committer=b"Test <test@example.com>",
    ).decode()

    (origin / "content.txt").write_text("version: beta\n")
    porcelain.add(str(origin), paths=["content.txt"])
    sha_beta = porcelain.commit(
        str(origin),
        message=b"beta commit",
        author=b"Test <test@example.com>",
        committer=b"Test <test@example.com>",
    ).decode()

    # Clone origin into the "cache" location (simulates get_or_update_cached_repo)
    cache_path = tmp_path / "cache" / "example.com" / "test" / "repo"
    cache_path.parent.mkdir(parents=True)
    porcelain.clone(str(origin), str(cache_path), checkout=True)

    return cache_path, sha_alpha, sha_beta


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestCheckoutNoCacheMutation:
    """
    checkout_to_work_dir must write files from the object store directly to
    work_dir and must never modify the cache repository's working tree or HEAD.
    """

    def test_correct_content_for_each_commit(self, two_commit_cache, tmp_path):
        """Sequential checkouts of two different commits yield the right files."""
        cache_path, sha_alpha, sha_beta = two_commit_cache
        fake_url = "https://example.com/test/repo"

        work_alpha = tmp_path / "work" / "alpha"
        work_beta = tmp_path / "work" / "beta"

        with patch(
            "omnibenchmark.git.cache.get_or_update_cached_repo",
            return_value=cache_path,
        ):
            _, resolved_alpha = checkout_to_work_dir(fake_url, sha_alpha, work_alpha)
            _, resolved_beta = checkout_to_work_dir(fake_url, sha_beta, work_beta)

        assert resolved_alpha == sha_alpha
        assert resolved_beta == sha_beta

        assert (work_alpha / "content.txt").read_text() == "version: alpha\n"
        assert (work_alpha / "extra.txt").read_text() == "extra file\n"
        assert (work_beta / "content.txt").read_text() == "version: beta\n"
        assert (work_beta / "extra.txt").read_text() == "extra file\n"

    def test_no_git_metadata_in_work_dir(self, two_commit_cache, tmp_path):
        """The work directory must not contain a .git directory."""
        cache_path, sha_alpha, _ = two_commit_cache
        fake_url = "https://example.com/test/repo"
        work_dir = tmp_path / "work" / "nogit"

        with patch(
            "omnibenchmark.git.cache.get_or_update_cached_repo",
            return_value=cache_path,
        ):
            checkout_to_work_dir(fake_url, sha_alpha, work_dir)

        assert not (work_dir / ".git").exists()

    def test_cache_head_not_changed_by_checkout(self, two_commit_cache, tmp_path):
        """
        Checking out the *earlier* commit must not reset the cache HEAD.

        Before the fix, checkout_to_work_dir called porcelain.reset on the
        shared cache repo, which mutated its working tree and HEAD.  After the
        fix, the cache must remain pointing at sha_beta (the latest commit).
        """
        cache_path, sha_alpha, sha_beta = two_commit_cache
        fake_url = "https://example.com/test/repo"

        cached_repo = Repo(str(cache_path))

        def _head_hex(repo: Repo) -> str:
            raw = repo.head()
            return raw.hex() if len(raw) == 20 else raw.decode("ascii")

        # Sanity: cache was cloned from origin whose HEAD is sha_beta
        head_before = _head_hex(cached_repo)
        assert head_before == sha_beta, "Fixture assumption violated"
        assert (cache_path / "content.txt").read_text() == "version: beta\n"

        # Check out the OLDER commit
        work_dir = tmp_path / "work" / "old"
        with patch(
            "omnibenchmark.git.cache.get_or_update_cached_repo",
            return_value=cache_path,
        ):
            checkout_to_work_dir(fake_url, sha_alpha, work_dir)

        # Work dir has the old content
        assert (work_dir / "content.txt").read_text() == "version: alpha\n"

        # Cache is completely unchanged
        head_after = _head_hex(cached_repo)
        assert (
            head_after == sha_beta
        ), f"Cache HEAD was modified!  Expected {sha_beta}, got {head_after}"
        assert (
            cache_path / "content.txt"
        ).read_text() == "version: beta\n", "Cache working tree was modified!"

    def test_short_hash_resolution(self, two_commit_cache, tmp_path):
        """A 7-character short hash is expanded to the full commit SHA."""
        cache_path, sha_alpha, _ = two_commit_cache
        fake_url = "https://example.com/test/repo"
        work_dir = tmp_path / "work" / "short"

        with patch(
            "omnibenchmark.git.cache.get_or_update_cached_repo",
            return_value=cache_path,
        ):
            _, resolved = checkout_to_work_dir(fake_url, sha_alpha[:7], work_dir)

        assert resolved == sha_alpha

    def test_unknown_ref_raises(self, two_commit_cache, tmp_path):
        """A ref that does not exist in the cache raises RuntimeError."""
        cache_path, _, _ = two_commit_cache
        fake_url = "https://example.com/test/repo"
        work_dir = tmp_path / "work" / "missing"

        with patch(
            "omnibenchmark.git.cache.get_or_update_cached_repo",
            return_value=cache_path,
        ):
            with pytest.raises(RuntimeError, match="not found in repository"):
                checkout_to_work_dir(fake_url, "deadbeef1234", work_dir)

    def test_overwrites_existing_work_dir(self, two_commit_cache, tmp_path):
        """A stale work_dir is replaced on re-checkout."""
        cache_path, sha_alpha, sha_beta = two_commit_cache
        fake_url = "https://example.com/test/repo"
        work_dir = tmp_path / "work" / "reuse"

        with patch(
            "omnibenchmark.git.cache.get_or_update_cached_repo",
            return_value=cache_path,
        ):
            checkout_to_work_dir(fake_url, sha_beta, work_dir)
            assert (work_dir / "content.txt").read_text() == "version: beta\n"

            # Re-checkout to earlier commit
            checkout_to_work_dir(fake_url, sha_alpha, work_dir)

        assert (work_dir / "content.txt").read_text() == "version: alpha\n"

    def test_concurrent_checkouts_of_different_refs(self, two_commit_cache, tmp_path):
        """
        N concurrent checkouts of two different commits from the same cache
        must all produce correct, uncorrupted content.

        This is the main regression test for the cache-mutation bug: the old
        code called porcelain.reset on the shared cache before copying, so
        a race between two threads could leave one work_dir with the wrong
        commit's files.
        """
        cache_path, sha_alpha, sha_beta = two_commit_cache
        fake_url = "https://example.com/test/repo"

        N = 6  # concurrent pairs
        results: dict = {}
        errors: list = []

        def do_checkout(sha: str, label: str) -> None:
            try:
                work_dir = tmp_path / "concurrent" / label
                checkout_to_work_dir(fake_url, sha, work_dir)
                results[label] = (work_dir / "content.txt").read_text()
            except Exception as exc:
                errors.append((label, repr(exc)))

        threads = [
            threading.Thread(target=do_checkout, args=(sha_alpha, f"alpha_{i}"))
            for i in range(N)
        ] + [
            threading.Thread(target=do_checkout, args=(sha_beta, f"beta_{i}"))
            for i in range(N)
        ]
        random.shuffle(threads)

        with patch(
            "omnibenchmark.git.cache.get_or_update_cached_repo",
            return_value=cache_path,
        ):
            for t in threads:
                t.start()
            for t in threads:
                t.join()

        assert not errors, f"Errors during concurrent checkout: {errors}"

        for i in range(N):
            assert (
                results.get(f"alpha_{i}") == "version: alpha\n"
            ), f"alpha_{i} had wrong content: {results.get(f'alpha_{i}')!r}"
            assert (
                results.get(f"beta_{i}") == "version: beta\n"
            ), f"beta_{i} had wrong content: {results.get(f'beta_{i}')!r}"
