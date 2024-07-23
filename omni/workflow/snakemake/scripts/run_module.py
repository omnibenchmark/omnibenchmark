#! /usr/bin/env python

import hashlib
import subprocess
import logging
import os
from typing import List

from git import Repo
from snakemake.script import Snakemake


def mock_execution(inputs: List[str], output: str, snakemake: Snakemake):
    print("Processed", inputs, "to", output, "using threads", snakemake.threads)
    print("  bench_iteration is", snakemake.bench_iteration)
    print("  resources are", snakemake.resources)
    print("  wildcards are", snakemake.wildcards)
    print("  rule is", snakemake.rule)
    print("  scriptdir is", snakemake.scriptdir)
    print("  params are", snakemake.params)


def execution(
    module_dir: str,
    module_name: str,
    output_dir: str,
    dataset: str,
    inputs_map: dict[str, str],
    parameters: List[str],
):
    run_sh = os.path.join(module_dir, "run.sh")
    if not os.path.exists(run_sh):
        logging.error(f"ERROR: {module_name} run.sh script does not exist.")
        raise RuntimeError(f"{module_name} run.sh script does not exist")

    # Constructing the command list
    command = [run_sh, output_dir, dataset]

    # Adding input files with their respective keys
    if inputs_map:
        for k, v in inputs_map.items():
            command.extend([f"--{k}", v])

    # Adding extra parameters
    if parameters:
        command.extend(parameters)

    try:
        # Execute the shell script
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        return result.stdout

    except subprocess.CalledProcessError as e:
        logging.error(
            f"ERROR: Executing {run_sh} failed with exit code {e.returncode}."
        )
        raise RuntimeError(
            f"ERROR: Executing {run_sh} failed with exit code {e.returncode}."
        ) from e


# Create a unique folder name based on the repository URL and commit hash
def generate_unique_repo_folder_name(repo_url, commit_hash):
    unique_string = f"{repo_url}@{commit_hash}"
    folder_name = hashlib.md5(unique_string.encode()).hexdigest()

    return folder_name


def clone_module(output_dir: str, repository_url: str, commit_hash: str):
    module_name = generate_unique_repo_folder_name(repository_url, commit_hash)
    module_dir = os.path.join(output_dir, module_name)

    if not os.path.exists(module_dir):
        logging.info(
            f"Cloning module `{repository_url}:{commit_hash}` to `{module_dir}`"
        )
        repo = Repo.clone_from(repository_url, module_dir)
        repo.git.checkout(commit_hash)
    else:
        repo = Repo(module_dir)

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
    repositories_dir = os.path.join(".snakemake", "repos")
    module_dir = clone_module(repositories_dir, repository_url, commit_hash)

    # Execute module code
    module_name = snakemake.rule

    output_dir = os.path.commonpath(snakemake.output)
    if len(snakemake.output) == 1:
        output_dir = os.path.dirname(output_dir)

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
