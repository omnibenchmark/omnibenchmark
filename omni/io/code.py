import hashlib
import logging
import os
from pathlib import Path
from typing import List

from git import Repo


# Create a unique folder name based on the repository URL and commit hash
def generate_unique_repo_folder_name(repo_url, commit_hash):
    unique_string = f"{repo_url}@{commit_hash}"
    folder_name = hashlib.md5(unique_string.encode()).hexdigest()

    return folder_name


def clone_module(output_dir: Path, repository_url: str, commit_hash: str) -> Path:
    module_name = generate_unique_repo_folder_name(repository_url, commit_hash)
    module_dir = output_dir / module_name

    if not module_dir.exists():
        logging.info(
            f"Cloning module `{repository_url}:{commit_hash}` to `{module_dir.as_posix()}`"
        )
        repo = Repo.clone_from(repository_url, module_dir.as_posix())
        repo.git.checkout(commit_hash)
    else:
        repo = Repo(module_dir.as_posix())

    if repo.head.commit.hexsha[:7] != commit_hash:
        logging.error(
            f"ERROR: Failed while cloning module `{repository_url}:{commit_hash}`"
        )
        logging.error(f"{commit_hash} does not match {repo.head.commit.hexsha[:7]}`")
        raise RuntimeError(
            f"ERROR: {commit_hash} does not match {repo.head.commit.hexsha[:7]}"
        )

    return module_dir


def dump_parameters_to_file(output_dir: str, parameters: List[str]):
    os.makedirs(output_dir, exist_ok=True)

    if parameters is not None:
        params_file = os.path.join(output_dir, "parameters.txt")
        with open(params_file, "w") as params_file:
            params_file.write(f"{parameters}")

        param_dict_file = os.path.join(output_dir, "..", "parameters_dict.txt")
        with open(param_dict_file, "a") as param_dict_file:
            param_dict_file.write(f"{os.path.basename(output_dir)} {parameters}\n")
