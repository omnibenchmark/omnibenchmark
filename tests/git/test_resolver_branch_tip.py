"""End-to-end regression for #186 at the resolver level.

`tests/git/test_cache_refresh.py` already proves the *cache* functions
(`checkout_to_work_dir`, `get_or_update_cached_repo`) follow a branch tip after a
push.  This module closes the remaining gap: it drives the full
``ModuleResolver.resolve`` path the way ``ob run`` does, against a **real bare
repo** updated with **`git push`**, and asserts that the *resolved commit*
(``ResolvedModule.commit``) and the checked-out entrypoint content advance to the
new tip on the second resolve — without clearing the cache.

Why the patches
---------------
``ModuleResolver`` treats any absolute filesystem path as a *local* repo
(``is_local_path``), which bypasses the git cache/fetch path entirely.  To
exercise the *remote* branch (cache → fetch → unpinned-ref resolution) we point
the module at a bare repo on disk but:

* patch ``resolver.is_local_path`` → ``False`` so the remote code path runs, and
* patch ``cache.parse_repo_url`` → a stable cache key (same trick as
  ``test_cache_refresh.py``).

Every other step — clone, fetch, ref resolution, work-dir checkout, entrypoint
dereference — runs unmodified.
"""

import shutil
import subprocess
from pathlib import Path

import pytest

from omnibenchmark.backend.resolver import ModuleResolver
from omnibenchmark.model.benchmark import Module, Repository

pytestmark = pytest.mark.skipif(
    shutil.which("git") is None, reason="git executable required"
)

_CACHE_KEY = "local/test/branch-tip-origin"
_ENV = {
    "GIT_AUTHOR_NAME": "Test",
    "GIT_AUTHOR_EMAIL": "test@example.com",
    "GIT_COMMITTER_NAME": "Test",
    "GIT_COMMITTER_EMAIL": "test@example.com",
}


def _git(cwd: Path, *args: str) -> str:
    """Run a git command, returning stripped stdout."""
    import os

    out = subprocess.run(
        ["git", *args],
        cwd=str(cwd),
        check=True,
        capture_output=True,
        text=True,
        env={**os.environ, **_ENV},
    )
    return out.stdout.strip()


def _write_module(work: Path, version: str) -> None:
    """Write a minimal host module whose entrypoint embeds *version*."""
    (work / "omnibenchmark.yaml").write_text("entrypoints:\n  default: run.sh\n")
    (work / "run.sh").write_text(
        "#!/usr/bin/env bash\n" f'VERSION="{version}"\n' 'echo "$VERSION"\n'
    )


@pytest.fixture()
def bare_remote(tmp_path):
    """A bare git repo (the "remote") with one pushed commit on ``main``.

    Returns ``(origin_bare, clone_dir, sha_v1)``.
    """
    origin = tmp_path / "origin.git"
    origin.mkdir()
    _git(origin, "init", "--bare", "--initial-branch=main")

    clone = tmp_path / "clone"
    clone.mkdir()
    _git(clone, "init", "--initial-branch=main")
    _git(clone, "remote", "add", "origin", str(origin))

    _write_module(clone, "v1")
    _git(clone, "add", "-A")
    _git(clone, "commit", "-m", "v1")
    _git(clone, "push", "-u", "origin", "main")
    sha_v1 = _git(clone, "rev-parse", "HEAD")

    return origin, clone, sha_v1


def _resolve(monkeypatch, resolver: ModuleResolver, module: Module):
    """Resolve a module through the *remote* path (patches applied)."""
    monkeypatch.setattr(
        "omnibenchmark.backend.resolver.is_local_path", lambda _url: False
    )
    monkeypatch.setattr(
        "omnibenchmark.git.cache.parse_repo_url", lambda _url: _CACHE_KEY
    )
    return resolver.resolve(
        module=module,
        module_id="M1",
        software_environment_id="host",
        unpinned=True,
    )


def _entrypoint_text(resolver: ModuleResolver, resolved) -> str:
    module_dir = Path(resolved.module_dir)
    if not module_dir.is_absolute():
        module_dir = resolver.output_dir / module_dir
    return (module_dir / resolved.entrypoint).read_text()


def test_resolver_picks_up_branch_tip_after_push(bare_remote, tmp_path, monkeypatch):
    """Second resolve of a branch ref must return the new commit and new content."""
    origin, clone, sha_v1 = bare_remote

    out_dir = tmp_path / "out"
    resolver = ModuleResolver(
        work_base_dir=out_dir / ".modules",
        cache_dir=tmp_path / "cache",
        output_dir=out_dir,
    )
    module = Module(
        id="M1",
        repository=Repository(url=str(origin), commit="main"),  # branch ref → #186
        software_environment="host",
    )

    # --- First resolve: sees v1 ---
    r1 = _resolve(monkeypatch, resolver, module)
    assert r1.commit == sha_v1
    assert 'VERSION="v1"' in _entrypoint_text(resolver, r1)

    # --- Push a new commit to the branch tip ---
    _write_module(clone, "v2")
    _git(clone, "commit", "-am", "v2")
    _git(clone, "push", "origin", "main")
    sha_v2 = _git(clone, "rev-parse", "HEAD")
    assert sha_v2 != sha_v1

    # --- Second resolve: must follow the tip (no cache clearing) ---
    r2 = _resolve(monkeypatch, resolver, module)
    assert r2.commit == sha_v2, "resolver served a stale cached commit (#186)"
    assert 'VERSION="v2"' in _entrypoint_text(resolver, r2)


def test_resolver_pinned_sha_unaffected_by_later_push(
    bare_remote, tmp_path, monkeypatch
):
    """Pinning the original SHA keeps resolving to v1 even after the branch moves."""
    origin, clone, sha_v1 = bare_remote

    out_dir = tmp_path / "out"
    resolver = ModuleResolver(
        work_base_dir=out_dir / ".modules",
        cache_dir=tmp_path / "cache",
        output_dir=out_dir,
    )

    # Advance the branch first so the cache must fetch newer objects.
    _write_module(clone, "v2")
    _git(clone, "commit", "-am", "v2")
    _git(clone, "push", "origin", "main")

    module = Module(
        id="M1",
        repository=Repository(url=str(origin), commit=sha_v1),  # pinned full SHA
        software_environment="host",
    )
    r = _resolve(monkeypatch, resolver, module)
    assert r.commit == sha_v1
    assert 'VERSION="v1"' in _entrypoint_text(resolver, r)
