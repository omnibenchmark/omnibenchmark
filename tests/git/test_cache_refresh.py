"""Integration tests for git cache refresh and branch-tracking behaviour.

These tests exercise the scenario that caused reported flakiness: a module
repository receives new commits after the local cache was first populated, and
the next resolution must return the *new* commit — not the stale cached one.

Strategy
--------
A local directory acts as the "remote" (dulwich can clone/fetch from plain
filesystem paths without a git server).  We patch only ``parse_repo_url`` so
that ``get_or_update_cached_repo`` maps our local path to a predictable cache
key; every other code path — the fetch, the object-store lookup, the
``build_index_from_tree`` write — runs unmodified.

Tests
-----
* ``test_checkout_reflects_new_commit_after_push``
      Main regression: second checkout returns the new SHA and new content.

* ``test_old_work_dir_unchanged_after_cache_update``
      The first work directory must not be mutated when the cache is refreshed.

* ``test_pinned_commit_still_resolves_after_cache_update``
      After new commits arrive, checking out the old pinned SHA still works and
      returns the original content.

* ``test_get_or_update_cached_repo_fetches_new_commits``
      Focused test of the fetch step: verifies that ``get_or_update_cached_repo``
      brings new objects into the cache so they can later be resolved.

* ``test_local_copy_tracks_head_after_new_commit``
      ``copy_local_to_work_dir`` (the non-cache, local-path code path) must also
      reflect HEAD after a new commit.
"""

import pytest
from unittest.mock import patch

from dulwich import porcelain
from dulwich.repo import Repo

from omnibenchmark.git.cache import (
    checkout_to_work_dir,
    copy_local_to_work_dir,
    get_or_update_cached_repo,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_AUTHOR = b"Test User <test@example.com>"
_CACHE_KEY = "local/test/origin"  # stable key returned by the patched parse_repo_url


def _commit(repo_dir, filename, content, message):
    """Write *filename* with *content*, stage it, and create a commit.

    Returns the 40-character hex SHA string of the new commit.
    """
    (repo_dir / filename).write_text(content)
    porcelain.add(str(repo_dir), paths=[filename])
    sha_bytes = porcelain.commit(
        str(repo_dir),
        message=message if isinstance(message, bytes) else message.encode(),
        author=_AUTHOR,
        committer=_AUTHOR,
    )
    return sha_bytes.decode("ascii")


def _head_sha(repo_dir):
    """Return the 40-char HEAD SHA of *repo_dir*."""
    with Repo(str(repo_dir)) as r:
        return r.head().decode("ascii")


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def local_remote_with_cache(tmp_path):
    """Return *(origin, cache_dir, fake_url, initial_sha)*.

    * ``origin``     — git repo used as the "remote"; commits here simulate
                       upstream pushes.
    * ``cache_dir``  — pre-populated cache (a clone of *origin* at
                       ``cache_dir / _CACHE_KEY``).
    * ``fake_url``   — the string passed as *repo_url* to the cache functions;
                       it is the local path of *origin* so dulwich can
                       clone/fetch from it, but ``parse_repo_url`` is patched
                       to map it to ``_CACHE_KEY``.
    * ``initial_sha``— 40-char hex SHA of the first commit in *origin*.
    """
    origin = tmp_path / "origin"
    origin.mkdir()
    porcelain.init(str(origin))
    # Ensure the default branch is named "main" regardless of system git config.
    with Repo(str(origin)) as _r:
        _r.refs.set_symbolic_ref(b"HEAD", b"refs/heads/main")
    initial_sha = _commit(origin, "sentinel.txt", "v1", b"initial: v1")

    # Populate cache exactly as get_or_update_cached_repo does on first call.
    cache_repo = tmp_path / "cache" / _CACHE_KEY
    cache_repo.parent.mkdir(parents=True)
    porcelain.clone(str(origin), str(cache_repo), checkout=True)

    return origin, tmp_path / "cache", str(origin), initial_sha


# ---------------------------------------------------------------------------
# Patch helper
# ---------------------------------------------------------------------------


def _patched_checkout(fake_url, ref, work_dir, cache_dir):
    """Call checkout_to_work_dir with parse_repo_url patched to _CACHE_KEY."""
    with patch(
        "omnibenchmark.git.cache.parse_repo_url",
        return_value=_CACHE_KEY,
    ):
        return checkout_to_work_dir(fake_url, ref, work_dir, cache_dir)


# ---------------------------------------------------------------------------
# Tests: cache refresh via two-tier system
# ---------------------------------------------------------------------------


class TestCacheRefreshTracksBranch:
    """checkout_to_work_dir must fetch new commits and reflect them in the work dir."""

    def test_checkout_reflects_new_commit_after_push(
        self, local_remote_with_cache, tmp_path
    ):
        """Second checkout on the same branch returns the new SHA and file content."""
        origin, cache_dir, fake_url, _initial_sha = local_remote_with_cache

        # First checkout — sees v1
        work_v1 = tmp_path / "work_v1"
        _, sha1 = _patched_checkout(fake_url, "main", work_v1, cache_dir)

        assert (work_v1 / "sentinel.txt").read_text() == "v1"

        # Push a new commit to origin (simulates upstream activity)
        sha2_expected = _commit(
            origin,
            "sentinel.txt",
            "SENTINEL=cache_refresh_test_v2",
            b"feat: sentinel v2",
        )

        # Second checkout — cache must fetch and resolve to the new commit
        work_v2 = tmp_path / "work_v2"
        _, sha2 = _patched_checkout(fake_url, "main", work_v2, cache_dir)

        assert sha2 != sha1
        assert sha2 == sha2_expected
        assert (
            work_v2 / "sentinel.txt"
        ).read_text() == "SENTINEL=cache_refresh_test_v2"

    def test_old_work_dir_unchanged_after_cache_update(
        self, local_remote_with_cache, tmp_path
    ):
        """The first work directory must not be mutated when the cache is refreshed."""
        origin, cache_dir, fake_url, _ = local_remote_with_cache

        work_v1 = tmp_path / "work_v1"
        _patched_checkout(fake_url, "main", work_v1, cache_dir)

        _commit(origin, "sentinel.txt", "v2 content", b"feat: v2")

        work_v2 = tmp_path / "work_v2"
        _patched_checkout(fake_url, "main", work_v2, cache_dir)

        # v1 work dir must still contain the original content
        assert (work_v1 / "sentinel.txt").read_text() == "v1"

    def test_pinned_commit_still_resolves_after_cache_update(
        self, local_remote_with_cache, tmp_path
    ):
        """After new commits arrive, the old pinned SHA still resolves correctly."""
        origin, cache_dir, fake_url, _ = local_remote_with_cache

        # Record sha1 via an initial checkout
        work_pin_ref = tmp_path / "work_pin_ref"
        _, sha1 = _patched_checkout(fake_url, "main", work_pin_ref, cache_dir)

        # Push more commits so the cache is refreshed on the next call
        _commit(origin, "sentinel.txt", "SENTINEL=newer_content", b"feat: newer")
        _commit(origin, "extra.txt", "extra", b"chore: add extra file")

        # Checkout the old pinned SHA — the old object must still be in the store
        work_pinned = tmp_path / "work_pinned"
        _, resolved = _patched_checkout(fake_url, sha1, work_pinned, cache_dir)

        assert resolved == sha1
        assert (work_pinned / "sentinel.txt").read_text() == "v1"
        assert not (work_pinned / "extra.txt").exists()


class TestGetOrUpdateCachedRepoFetches:
    """get_or_update_cached_repo must bring new objects into the cache."""

    def test_fetch_updates_remote_tracking_ref(self, local_remote_with_cache, tmp_path):
        """After a new commit on origin, the remote-tracking ref in the cache advances."""
        origin, cache_dir, fake_url, _ = local_remote_with_cache

        sha_new = _commit(origin, "sentinel.txt", "v2", b"feat: v2")

        with patch(
            "omnibenchmark.git.cache.parse_repo_url",
            return_value=_CACHE_KEY,
        ):
            cache_path = get_or_update_cached_repo(fake_url, cache_dir)

        # The remote-tracking branch in the cache should now point to sha_new
        with Repo(str(cache_path)) as cached_repo:
            remote_ref = b"refs/remotes/origin/main"
            assert remote_ref in cached_repo.refs
            cached_sha = cached_repo.refs[remote_ref].decode("ascii")

        assert cached_sha == sha_new


# ---------------------------------------------------------------------------
# Tests: local copy path (no cache involved)
# ---------------------------------------------------------------------------


class TestLocalCopyTracksHead:
    """copy_local_to_work_dir must always reflect the current HEAD of the local repo."""

    def test_reflects_new_commit_after_local_commit(self, tmp_path):
        """New commit in the local repo is picked up by a subsequent copy."""
        repo_dir = tmp_path / "local_repo"
        repo_dir.mkdir()
        porcelain.init(str(repo_dir))
        sha1 = _commit(repo_dir, "sentinel.txt", "v1", b"initial")

        work_v1 = tmp_path / "work_v1"
        _, resolved1 = copy_local_to_work_dir(repo_dir, work_v1)
        assert resolved1 == sha1
        assert (work_v1 / "sentinel.txt").read_text() == "v1"

        sha2 = _commit(
            repo_dir, "sentinel.txt", "SENTINEL=local_head_test_v2", b"feat: v2"
        )

        work_v2 = tmp_path / "work_v2"
        _, resolved2 = copy_local_to_work_dir(repo_dir, work_v2)

        assert resolved2 == sha2
        assert resolved2 != resolved1
        assert (work_v2 / "sentinel.txt").read_text() == "SENTINEL=local_head_test_v2"
        # v1 work dir is independent
        assert (work_v1 / "sentinel.txt").read_text() == "v1"
