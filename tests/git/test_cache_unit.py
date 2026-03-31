"""Short unit tests for git/cache.py — no real git operations, all dulwich mocked."""

import pytest
from pathlib import Path
from unittest.mock import MagicMock, patch

from dulwich.refs import Ref

from omnibenchmark.git.cache import checkout_to_work_dir, get_or_update_cached_repo


# ── helpers ────────────────────────────────────────────────────────────────────


def _commit_obj(hex_sha: str):
    """Minimal mock dulwich commit object whose .id is 20-byte binary."""
    obj = MagicMock()
    obj.id = bytes.fromhex(hex_sha)
    obj.tree = b"\x00" * 20
    return obj


_SHA1 = "a" * 40


def _checkout(ref, refs_dict, object_store_ids, parse_side_effect, tmp_path):
    """
    Run checkout_to_work_dir with all dulwich calls replaced by mocks.

    parse_side_effect is passed directly as the side_effect of the
    dulwich.objectspec.parse_commit mock, so it can be a list of
    return-values / exceptions or a single exception class.
    """
    mock_repo = MagicMock()
    mock_repo.refs = refs_dict
    mock_repo.object_store = object_store_ids

    work_dir = tmp_path / "work"
    cache_dir = tmp_path / "cache"

    with (
        patch(
            "omnibenchmark.git.cache.get_or_update_cached_repo",
            return_value=cache_dir,
        ),
        patch(
            "omnibenchmark.git.cache.porcelain.open_repo",
            return_value=mock_repo,
        ),
        patch(
            "dulwich.objectspec.parse_commit",
            side_effect=parse_side_effect,
        ),
        patch("dulwich.index.build_index_from_tree"),
    ):
        return checkout_to_work_dir("fake://url", ref, work_dir, cache_dir)


# ── ref resolution ─────────────────────────────────────────────────────────────


@pytest.mark.short
class TestCheckoutRefResolution:
    """Unit tests for the three-stage ref resolution in checkout_to_work_dir."""

    def test_remote_tracking_ref_resolves(self, tmp_path):
        """refs/remotes/origin/<branch> is used when present."""
        sha_bytes = bytes.fromhex(_SHA1)
        commit = _commit_obj(_SHA1)
        refs = {Ref(b"refs/remotes/origin/main"): sha_bytes}

        _, sha = _checkout("main", refs, [], [commit], tmp_path)

        assert sha == _SHA1

    def test_direct_parse_used_when_no_remote_ref(self, tmp_path):
        """Full SHA with no matching remote tracking ref resolves via direct parse."""
        commit = _commit_obj(_SHA1)

        _, sha = _checkout(_SHA1, {}, [], [commit], tmp_path)

        assert sha == _SHA1

    def test_direct_parse_keyerror_falls_through(self, tmp_path):
        """KeyError from direct parse is silently swallowed; short hash path runs next."""
        sha_bytes = bytes.fromhex(_SHA1)
        commit = _commit_obj(_SHA1)
        # First call (direct parse) raises KeyError; second call (short-hash) succeeds.
        _, sha = _checkout(
            _SHA1[:8],
            {},
            [sha_bytes],
            [KeyError, commit],
            tmp_path,
        )

        assert sha == _SHA1

    def test_short_hash_expansion(self, tmp_path):
        """A 7-char hex prefix resolves via object_store scan."""
        sha_bytes = bytes.fromhex(_SHA1)
        commit = _commit_obj(_SHA1)

        _, sha = _checkout(
            _SHA1[:7],
            {},
            [sha_bytes],
            [KeyError, commit],
            tmp_path,
        )

        assert sha == _SHA1

    def test_ambiguous_short_hash_raises(self, tmp_path):
        """Two objects sharing the same prefix → RuntimeError 'ambiguous'."""
        # Both SHAs start with "abcd" (≥4 chars required by the regex guard)
        sha_a = "abcd" + "00" * 18
        sha_b = "abcd" + "ff" * 18

        with pytest.raises(RuntimeError, match="ambiguous"):
            _checkout(
                "abcd",
                {},
                [bytes.fromhex(sha_a), bytes.fromhex(sha_b)],
                KeyError,
                tmp_path,
            )

    def test_unknown_ref_raises(self, tmp_path):
        """A ref that matches nothing raises RuntimeError 'not found'."""
        with pytest.raises(RuntimeError, match="not found"):
            _checkout(
                "nonexistent-branch",
                {},
                [],
                KeyError,
                tmp_path,
            )

    def test_existing_work_dir_is_replaced(self, tmp_path):
        """If work_dir already exists its contents are replaced by the checkout."""
        work_dir = tmp_path / "work"
        work_dir.mkdir()
        stale = work_dir / "stale.txt"
        stale.write_text("old")

        sha_bytes = bytes.fromhex(_SHA1)
        commit = _commit_obj(_SHA1)
        refs = {Ref(b"refs/remotes/origin/main"): sha_bytes}

        mock_repo = MagicMock()
        mock_repo.refs = refs
        mock_repo.object_store = []
        cache_dir = tmp_path / "cache"

        with (
            patch(
                "omnibenchmark.git.cache.get_or_update_cached_repo",
                return_value=cache_dir,
            ),
            patch(
                "omnibenchmark.git.cache.porcelain.open_repo",
                return_value=mock_repo,
            ),
            patch("dulwich.objectspec.parse_commit", return_value=commit),
            patch("dulwich.index.build_index_from_tree"),
        ):
            checkout_to_work_dir("fake://url", "main", work_dir, cache_dir)

        assert not stale.exists()


# ── corrupt cache fallback ─────────────────────────────────────────────────────


@pytest.mark.short
class TestGetOrUpdateCachedRepo:
    """Unit tests for get_or_update_cached_repo error paths."""

    def test_corrupt_cache_triggers_reclone(self, tmp_path):
        """If fetch raises, the corrupt dir is removed and the repo is re-cloned."""
        cache_dir = tmp_path / "cache"
        repo_cache_dir = cache_dir / "local" / "test" / "origin"
        repo_cache_dir.mkdir(parents=True)
        sentinel = repo_cache_dir / "sentinel"
        sentinel.write_text("corrupt")

        with (
            patch(
                "omnibenchmark.git.cache.parse_repo_url",
                return_value="local/test/origin",
            ),
            patch(
                "omnibenchmark.git.cache.porcelain.open_repo",
                side_effect=Exception("corrupt repo"),
            ),
            patch("omnibenchmark.git.cache.porcelain.clone") as mock_clone,
        ):
            get_or_update_cached_repo("fake://url", cache_dir)

        # The corrupt directory must have been wiped before re-cloning.
        assert not sentinel.exists()
        mock_clone.assert_called_once()

    def test_clone_failure_cleans_partial_dir_and_raises(self, tmp_path):
        """If the clone partially creates the target dir then fails, it is cleaned up."""
        cache_dir = tmp_path / "cache"
        target = cache_dir / "local" / "test" / "origin"

        def _partial_clone(source, target, **kwargs):
            Path(target).mkdir(parents=True, exist_ok=True)
            (Path(target) / "partial").write_text("debris")
            raise Exception("network error")

        with (
            patch(
                "omnibenchmark.git.cache.parse_repo_url",
                return_value="local/test/origin",
            ),
            patch(
                "omnibenchmark.git.cache.porcelain.clone",
                side_effect=_partial_clone,
            ),
        ):
            with pytest.raises(RuntimeError, match="Failed to clone"):
                get_or_update_cached_repo("fake://url", cache_dir)

        assert not target.exists()
