import logging
import re
import subprocess
from typing import Optional

from pathlib import Path
from filelock import FileLock
from dulwich import porcelain
from dulwich.errors import NotGitRepository

from omnibenchmark.config import get_git_modules_dir
from omnibenchmark.model.repo import get_repo_hash


def is_commit_hash(ref: str) -> bool:
    """
    Check if a reference looks like a commit hash.

    Args:
        ref: Git reference (could be branch, tag, or commit hash)

    Returns:
        True if ref appears to be a commit hash, False otherwise
    """
    # Commit hashes are typically 7-40 hex characters
    return bool(re.match(r"^[a-f0-9]{7,40}$", ref, re.IGNORECASE))


def _get_head_sha(repo_path: Path) -> str:
    """Return the short (7-char) HEAD commit SHA for a local repo."""
    result = subprocess.run(
        ["git", "-C", str(repo_path), "rev-parse", "--short", "HEAD"],
        capture_output=True,
        text=True,
        check=True,
    )
    return result.stdout.strip()


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
    try:
        # If it looks like a branch name, try to pull latest changes
        if not is_commit_hash(commit_or_branch):
            logging.info(
                f"Detected branch name '{commit_or_branch}', pulling latest changes"
            )
            try:
                subprocess.run(
                    ["git", "-C", str(repo_path), "fetch", "origin"],
                    capture_output=True,
                    check=True,
                )
                subprocess.run(
                    ["git", "-C", str(repo_path), "checkout", commit_or_branch],
                    capture_output=True,
                    check=True,
                )
                subprocess.run(
                    ["git", "-C", str(repo_path), "pull", "origin", commit_or_branch],
                    capture_output=True,
                    check=True,
                )
                logging.info(
                    f"Updated local branch '{commit_or_branch}' with latest changes"
                )
            except subprocess.CalledProcessError as e:
                logging.warning(
                    f"Failed to pull branch '{commit_or_branch}': {e}. Treating as commit hash."
                )
                subprocess.run(
                    ["git", "-C", str(repo_path), "checkout", commit_or_branch],
                    capture_output=True,
                    check=True,
                )
        else:
            # It's a commit hash, just checkout
            subprocess.run(
                ["git", "-C", str(repo_path), "checkout", commit_or_branch],
                capture_output=True,
                check=True,
            )

        return _get_head_sha(repo_path)

    except Exception as e:
        logging.error(f"Failed to update repository: {e}")
        # Return current commit as fallback
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
    porcelain.clone(repository_url, str(output_dir), checkout=False)
    subprocess.run(
        ["git", "-C", str(output_dir), "checkout", commit_hash],
        capture_output=True,
        check=True,
    )
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
                porcelain.clone(repository_url, str(module_dir))
                subprocess.run(
                    ["git", "-C", str(module_dir), "checkout", commit_or_branch],
                    capture_output=True,
                    check=True,
                )
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
