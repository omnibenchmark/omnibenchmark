#! /usr/bin/env python

import hashlib
import logging
import os
from pathlib import Path
from typing import List

from git import Repo
from snakemake.script import Snakemake

from omni.workflow.snakemake.scripts.execution import execution


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


try:
    snakemake: Snakemake = snakemake
    params = dict(snakemake.params)

    parameters = params["parameters"]
    repository_url = params["repository_url"]
    commit_hash = params["commit_hash"]
    inputs_map = params.get("inputs_map")
    dataset = params.get("dataset")
    if dataset is None:
        dataset = getattr(snakemake.wildcards, "dataset", "unknown")

    # Create parameters file for outputs
    output_dir = os.path.dirname(snakemake.output[0])
    dump_parameters_to_file(output_dir, parameters)

    # Clone github repository
    repositories_dir = Path(".snakemake") / "repos"
    module_dir = clone_module(repositories_dir, repository_url, commit_hash)

    # Execute module code
    module_name = snakemake.rule

    output_dir = Path(str(os.path.commonpath(snakemake.output)))
    if len(snakemake.output) == 1:
        output_dir = Path(os.path.dirname(output_dir))

    execution(
        module_dir,
        module_name=module_name,
        output_dir=output_dir,
        dataset=dataset,
        inputs_map=inputs_map,
        parameters=parameters,
    )

except NameError:
    raise RuntimeError(f"This script must be run from within a Snakemake workflow")
