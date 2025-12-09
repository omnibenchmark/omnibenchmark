import logging

from pathlib import Path
from filelock import FileLock
from git.exc import InvalidGitRepositoryError
from git import Repo

from omnibenchmark.workflow.snakemake.scripts.utils import (
    generate_unique_repo_folder_name,
)


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


def clone_module(output_dir: Path, repository_url: str, commit_hash: str) -> Path:
    module_name = generate_unique_repo_folder_name(repository_url, commit_hash)
    module_dir = output_dir / module_name

    lock = module_dir.with_suffix(".lock")
    with FileLock(lock):
        if not module_dir.exists():
            try:
                logging.info(
                    f"Get git archive `{repository_url}:{commit_hash}` to `{module_dir.as_posix()}`"
                )
                clone_git_repo(module_dir, repository_url, commit_hash)
                observed_commit_hash = commit_hash
            except Exception:
                logging.info(
                    f"Archival retrieval failed, cloning module `{repository_url}:{commit_hash}` to `{module_dir.as_posix()}`"
                )
                repo = Repo.clone_from(repository_url, module_dir.as_posix())
                repo.git.checkout(commit_hash)
                observed_commit_hash = repo.head.commit.hexsha[:7]
        else:
            try:
                repo = Repo(module_dir.as_posix())
                observed_commit_hash = repo.head.commit.hexsha[:7]
            except InvalidGitRepositoryError:
                # is archive
                observed_commit_hash = commit_hash
            except Exception:
                observed_commit_hash = "known"

        if observed_commit_hash != commit_hash:
            logging.error(
                f"ERROR: Failed while cloning module `{repository_url}:{commit_hash}`"
            )
            logging.error(
                f"{commit_hash} does not match {repo.head.commit.hexsha[:7]}`"
            )
            raise RuntimeError(
                f"ERROR: {commit_hash} does not match {repo.head.commit.hexsha[:7]}"
            )

    return module_dir
