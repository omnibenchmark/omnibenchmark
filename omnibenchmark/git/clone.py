import io
import logging
import re
from typing import Optional, cast
from pathlib import Path

from filelock import FileLock
from dulwich import porcelain
from dulwich.errors import NotGitRepository
from dulwich.objectspec import parse_commit
from dulwich.refs import Ref
from dulwich.repo import Repo

from omnibenchmark.config import get_git_modules_dir
from omnibenchmark.model.repo import get_repo_hash


class _NullBinaryStream(io.BytesIO):
    def write(self, b) -> int:
        return len(b) if b else 0


class _NullTextStream(io.StringIO):
    def write(self, s: str) -> int:
        return len(s)


_DEVNULL = _NullBinaryStream()
_DEVNULL_TEXT = _NullTextStream()


def is_commit_hash(ref: str) -> bool:
    """
    Check if a reference looks like a commit hash.

    Args:
        ref: Git reference (could be branch, tag, or commit hash)

    Returns:
        True if ref appears to be a commit hash, False otherwise
    """
    return bool(re.match(r"^[a-f0-9]{7,40}$", ref, re.IGNORECASE))


def _get_head_sha(repo_path: Path) -> str:
    """Return the short (7-char) HEAD commit SHA for a local repo."""
    repo = cast(Repo, porcelain.open_repo(str(repo_path)))
    head_bytes = repo.head()
    full_hex = head_bytes.hex() if len(head_bytes) == 20 else head_bytes.decode("ascii")
    return full_hex[:7]


def _check_clean_or_raise(repo_path: Path) -> None:
    """
    Raise a fatal error if the repository has local modifications.

    The module cache must never be edited by hand. If staged or unstaged
    changes are detected, we refuse to continue rather than silently
    overwriting the user's work.
    """
    repo = cast(Repo, porcelain.open_repo(str(repo_path)))
    status = porcelain.status(repo)

    dirty: list[bytes] = []
    for paths in status.staged.values():
        dirty.extend(paths)
    dirty.extend(status.unstaged)

    if dirty:
        files = "\n  ".join(
            p.decode() if isinstance(p, bytes) else p for p in dirty[:20]
        )
        raise RuntimeError(
            f"FATAL: module cache at '{repo_path}' has local modifications.\n"
            f"Modified files:\n  {files}\n\n"
            "The module cache must not be edited manually. "
            "Remove the directory and let omnibenchmark re-clone it:\n"
            f"  rm -rf '{repo_path}'"
        )


def _resolve_and_reset(repo: Repo, repo_path: Path, commit_or_branch: str) -> None:
    """
    Resolve commit_or_branch to a commit object and hard-reset the working tree.
    Mirrors the ref resolution logic from cache.py.
    """
    commit_obj = None

    # 1. Try direct resolution (full hash, tag, local branch)
    try:
        commit_obj = parse_commit(repo, commit_or_branch.encode("ascii"))
    except KeyError:
        pass

    # 2. Try remote tracking branch: refs/remotes/origin/<branch>
    if commit_obj is None:
        remote_ref = Ref(f"refs/remotes/origin/{commit_or_branch}".encode("ascii"))
        if remote_ref in repo.refs:
            sha = repo.refs[remote_ref]
            commit_obj = parse_commit(repo, sha)

    # 3. Expand short hashes
    if commit_obj is None and re.match(
        r"^[0-9a-f]{4,39}$", commit_or_branch, re.IGNORECASE
    ):
        prefix = commit_or_branch.lower()
        matches = [
            obj_id.hex()
            for obj_id in repo.object_store
            if obj_id.hex().startswith(prefix)
        ]
        if len(matches) == 1:
            commit_obj = parse_commit(repo, matches[0].encode("ascii"))
        elif len(matches) > 1:
            raise RuntimeError(
                f"Reference '{commit_or_branch}' is ambiguous "
                f"(matches {len(matches)} objects) in '{repo_path}'"
            )

    if commit_obj is None:
        raise RuntimeError(
            f"Reference '{commit_or_branch}' not found in repository at '{repo_path}'"
        )

    porcelain.reset(repo, "hard", commit_obj)


def update_existing_repo(
    repo_path: Path, repository_url: str, commit_or_branch: str
) -> str:
    """
    Update an existing repository, pulling latest changes if dealing with a branch.

    Args:
        repo_path: Path to the local git repository
        repository_url: URL of the git repository
        commit_or_branch: Either a commit hash or branch name

    Returns:
        The actual commit hash after update
    """
    _check_clean_or_raise(repo_path)

    try:
        repo = cast(Repo, porcelain.open_repo(str(repo_path)))

        if not is_commit_hash(commit_or_branch):
            logging.info(
                f"Detected branch name '{commit_or_branch}', pulling latest changes"
            )
            porcelain.fetch(
                repo, repository_url, errstream=_DEVNULL, outstream=_DEVNULL_TEXT
            )
            logging.info(
                f"Updated local branch '{commit_or_branch}' with latest changes"
            )

        _resolve_and_reset(repo, repo_path, commit_or_branch)
        return _get_head_sha(repo_path)

    except RuntimeError:
        raise
    except Exception as e:
        logging.error(f"Failed to update repository: {e}")
        try:
            return _get_head_sha(repo_path)
        except Exception:
            return "unknown"


def clone_git_repo(output_dir: Path, repository_url: str, commit_hash: str) -> Path:
    """
    Perform a shallow clone of a git repository at a specific commit.

    Args:
        output_dir: Directory where to clone the repository
        repository_url: URL of the git repository
        commit_hash: Specific commit hash to checkout
    """
    porcelain.clone(
        repository_url,
        str(output_dir),
        checkout=False,
        errstream=_DEVNULL,
        outstream=_DEVNULL,
    )
    repo = cast(Repo, porcelain.open_repo(str(output_dir)))
    _resolve_and_reset(repo, output_dir, commit_hash)
    return output_dir


def clone_module(
    repository_url: str, commit_or_branch: str, output_dir: Optional[Path] = None
) -> Path:
    """
    Clone or update a git module, handling both commit hashes and branch names.

    For branch names, will pull latest changes if the repository already exists.
    For commit hashes, will checkout the specific commit.

    Args:
        repository_url: URL of the git repository
        commit_or_branch: Either a specific commit hash or branch name
        output_dir: Directory where to clone the repository (defaults to configured git modules dir)

    Returns:
        Path to the cloned/updated module directory
    """
    if output_dir is None:
        output_dir = get_git_modules_dir()

    module_name = get_repo_hash(repository_url, commit_or_branch)
    module_dir = output_dir / module_name

    lock = module_dir.with_suffix(".lock")
    with FileLock(lock):
        if not module_dir.exists():
            try:
                logging.info(
                    f"Get git archive `{repository_url}:{commit_or_branch}` to `{module_dir.as_posix()}`"
                )
                clone_git_repo(module_dir, repository_url, commit_or_branch)
                observed_commit_hash = (
                    commit_or_branch if is_commit_hash(commit_or_branch) else None
                )
            except Exception:
                logging.info(
                    f"Archival retrieval failed, cloning module `{repository_url}:{commit_or_branch}` to `{module_dir.as_posix()}`"
                )
                porcelain.clone(
                    repository_url,
                    str(module_dir),
                    errstream=_DEVNULL,
                    outstream=_DEVNULL,
                )
                repo = cast(Repo, porcelain.open_repo(str(module_dir)))
                _resolve_and_reset(repo, module_dir, commit_or_branch)
                observed_commit_hash = _get_head_sha(module_dir)
        else:
            try:
                # Verify it's a git repo
                porcelain.open_repo(str(module_dir))
                # Update repository if it exists, handling branches vs commits
                observed_commit_hash = update_existing_repo(
                    module_dir, repository_url, commit_or_branch
                )
            except NotGitRepository:
                # is archive - no git repo, can't update
                observed_commit_hash = commit_or_branch
                logging.info(
                    f"Directory {module_dir} exists but is not a git repository (likely an archive)"
                )
            except RuntimeError:
                raise
            except Exception as e:
                logging.warning(f"Failed to update existing repository: {e}")
                observed_commit_hash = "known"

        # Only validate commit hash match for actual commit hashes, not branch names
        if (
            is_commit_hash(commit_or_branch)
            and observed_commit_hash != commit_or_branch
        ):
            logging.error(
                f"ERROR: Failed while cloning module `{repository_url}:{commit_or_branch}`"
            )
            logging.error(f"{commit_or_branch} does not match {observed_commit_hash}")
            raise RuntimeError(
                f"ERROR: {commit_or_branch} does not match {observed_commit_hash}"
            )
        elif not is_commit_hash(commit_or_branch):
            logging.info(
                f"Successfully updated branch '{commit_or_branch}' to commit {observed_commit_hash}"
            )

    return module_dir
