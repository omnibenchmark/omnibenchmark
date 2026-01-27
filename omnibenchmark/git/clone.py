import logging
import re
from typing import Optional

from pathlib import Path
from filelock import FileLock
from git.exc import InvalidGitRepositoryError, GitCommandError
from git import Repo

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


def update_existing_repo(repo: Repo, repository_url: str, commit_or_branch: str) -> str:
    """
    Update an existing repository, pulling latest changes if dealing with a branch.

    Args:
        repo: Git repository object
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
                # Fetch latest changes from remote
                repo.remotes.origin.fetch()

                # Check if the branch exists locally
                try:
                    local_branch = repo.heads[commit_or_branch]
                    # Switch to the branch and pull
                    local_branch.checkout()
                    repo.remotes.origin.pull(commit_or_branch)
                    logging.info(
                        f"Updated local branch '{commit_or_branch}' with latest changes"
                    )
                except IndexError:
                    # Branch doesn't exist locally, create it tracking remote
                    remote_branch = f"origin/{commit_or_branch}"
                    if remote_branch in [ref.name for ref in repo.remotes.origin.refs]:
                        repo.create_head(commit_or_branch, f"origin/{commit_or_branch}")
                        repo.heads[commit_or_branch].set_tracking_branch(
                            repo.remotes.origin.refs[commit_or_branch]
                        )
                        repo.heads[commit_or_branch].checkout()
                        logging.info(
                            f"Created and checked out tracking branch '{commit_or_branch}'"
                        )
                    else:
                        logging.warning(
                            f"Remote branch '{commit_or_branch}' not found, treating as commit hash"
                        )
                        repo.git.checkout(commit_or_branch)

            except GitCommandError as e:
                logging.warning(
                    f"Failed to pull branch '{commit_or_branch}': {e}. Treating as commit hash."
                )
                repo.git.checkout(commit_or_branch)
        else:
            # It's a commit hash, just checkout
            repo.git.checkout(commit_or_branch)

        return repo.head.commit.hexsha[:7]

    except Exception as e:
        logging.error(f"Failed to update repository: {e}")
        # Return current commit as fallback
        return repo.head.commit.hexsha[:7]


def clone_git_repo(output_dir: Path, repository_url: str, commit_hash: str) -> Path:
    """
    Perform a shallow clone of a git repository at a specific commit.

    Args:
        output_dir: Directory where to clone the repository
        repository_url: URL of the git repository
        commit_hash: Specific commit hash to checkout
    """
    repo = Repo.clone_from(repository_url, output_dir.as_posix(), no_checkout=True)
    repo.git.checkout(commit_hash)
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
                repo = Repo.clone_from(repository_url, module_dir.as_posix())
                repo.git.checkout(commit_or_branch)
                observed_commit_hash = repo.head.commit.hexsha[:7]
        else:
            try:
                repo = Repo(module_dir.as_posix())
                # Update repository if it exists, handling branches vs commits
                observed_commit_hash = update_existing_repo(
                    repo, repository_url, commit_or_branch
                )
            except InvalidGitRepositoryError:
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
