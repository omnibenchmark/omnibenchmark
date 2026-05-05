"""Regression tests for ModuleResolver concurrency.

The resolver runs under ThreadPoolExecutor(max_workers=cores) in cli/run.py.
When multiple modules share the same local repository (e.g. two modules
with `url: .` in the same benchmark), the local-path branch of
``_resolve_repository`` previously raced on the rmtree+copytree+move of
the shared ``pending/`` working directory, producing intermittent
``[Errno 17] File exists`` or ``[Errno 2] No such file`` failures during
``ob run --dirty --cores N``.

These tests pin the fix: a per-resolver lock around the local-path
read-decide-copy-write block so the existing ``_dirty_copy_cache`` can
deduplicate copies as intended.
"""

from __future__ import annotations

import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pytest

from omnibenchmark.backend.resolver import ModuleResolver
from omnibenchmark.model import Module, Repository


def _init_local_module_repo(repo_dir: Path) -> None:
    """Create a minimal git repo with an omnibenchmark.yaml + entrypoint."""
    repo_dir.mkdir(parents=True, exist_ok=True)
    (repo_dir / "omnibenchmark.yaml").write_text("entrypoints:\n  default: run.sh\n")
    (repo_dir / "run.sh").write_text("#!/usr/bin/env bash\necho hi\n")
    (repo_dir / "run.sh").chmod(0o755)

    def git(*args: str) -> None:
        subprocess.run(
            ["git", *args],
            cwd=repo_dir,
            check=True,
            capture_output=True,
        )

    git("init", "-q", "-b", "main")
    git("config", "user.email", "test@example.com")
    git("config", "user.name", "test")
    git("add", ".")
    git("commit", "-q", "-m", "init")


@pytest.fixture
def local_module_repo(tmp_path: Path) -> Path:
    repo_dir = tmp_path / "myrepo"
    _init_local_module_repo(repo_dir)
    return repo_dir


def _make_module(local_path: Path, module_id: str) -> Module:
    return Module(
        id=module_id,
        repository=Repository(
            url=str(local_path),
            commit="main",
            entrypoint="default",
        ),
        software_environment="env",
    )


@pytest.mark.short
class TestConcurrentLocalPathResolve:
    """Two modules sharing a local repo must resolve concurrently without racing."""

    def test_resolver_has_dirty_copy_lock(self, tmp_path: Path) -> None:
        """The fix exposes a lock on the resolver — sanity check it exists."""
        resolver = ModuleResolver(work_base_dir=tmp_path / "work")
        assert hasattr(resolver, "_dirty_copy_lock"), (
            "ModuleResolver must hold a lock to serialize local-path "
            "rmtree/copytree under ThreadPoolExecutor"
        )

    def test_concurrent_resolve_same_local_repo(
        self, tmp_path: Path, local_module_repo: Path
    ) -> None:
        """Two threads resolving modules sharing one local repo must both succeed."""
        resolver = ModuleResolver(
            work_base_dir=tmp_path / "work",
            output_dir=tmp_path / "out",
        )
        modules = [
            _make_module(local_module_repo, "module_a"),
            _make_module(local_module_repo, "module_b"),
        ]

        with ThreadPoolExecutor(max_workers=2) as pool:
            futures = [
                pool.submit(
                    resolver.resolve,
                    module=m,
                    module_id=m.id,
                    software_environment_id=m.software_environment,
                    dirty=True,
                )
                for m in modules
            ]
            resolved = [f.result() for f in as_completed(futures)]

        # Both resolutions should land on the same shared work dir — the
        # cache deduplicates copies of the same local repo within one
        # resolver lifetime.
        work_dirs = {r.module_dir for r in resolved}
        assert (
            len(work_dirs) == 1
        ), f"Expected both modules to share one work dir, got {work_dirs}"

        # Cache must be populated with exactly one entry keyed on the
        # local repo path.
        assert list(resolver._dirty_copy_cache.keys()) == [str(local_module_repo)]

    def test_high_concurrency_resolve_same_local_repo(
        self, tmp_path: Path, local_module_repo: Path
    ) -> None:
        """Stress: many parallel resolves of one local repo, no exceptions."""
        resolver = ModuleResolver(
            work_base_dir=tmp_path / "work",
            output_dir=tmp_path / "out",
        )
        n = 8
        modules = [_make_module(local_module_repo, f"module_{i}") for i in range(n)]

        with ThreadPoolExecutor(max_workers=n) as pool:
            futures = [
                pool.submit(
                    resolver.resolve,
                    module=m,
                    module_id=m.id,
                    software_environment_id=m.software_environment,
                    dirty=True,
                )
                for m in modules
            ]
            resolved = [f.result() for f in as_completed(futures)]

        assert len(resolved) == n
        work_dirs = {r.module_dir for r in resolved}
        assert len(work_dirs) == 1
