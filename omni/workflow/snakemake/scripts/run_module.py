#! /usr/bin/env python
# WARNING: Custom dependencies might not be available here, since this is run inside a specified environment.

import hashlib
import logging
import os
import subprocess
import sys
import csv
from datetime import datetime
from pathlib import Path
from typing import List

from snakemake.script import Snakemake

from omni.workflow.snakemake.scripts.execution import execution
from omni.workflow.snakemake.scripts.utils import *


# Create a unique folder name based on the repository URL and commit hash
def generate_unique_repo_folder_name(repo_url: str, commit_hash: str) -> str:
    unique_string = f"{repo_url}@{commit_hash}"
    folder_name = hashlib.md5(unique_string.encode()).hexdigest()

    return folder_name


def generate_metric_id(metric_name: str, output_path: Path, output_file: str) -> str:
    unique_string = f"{metric_name}@{output_path}@{output_file}"
    metric_id = hashlib.md5(unique_string.encode()).hexdigest()

    return metric_id


def clone_module(output_dir: Path, repository_url: str, commit_hash: str) -> Path:
    try:
        import git  # Check if gitpython is available
    except ImportError:
        print("gitpython is not installed. Installing now...")
        subprocess.check_call(
            [sys.executable, "-m", "pip", "install", "gitpython==3.1.43"]
        )
        import git  # Retry importing

    module_name = generate_unique_repo_folder_name(repository_url, commit_hash)
    module_dir = output_dir / module_name

    if not module_dir.exists():
        logging.info(
            f"Cloning module `{repository_url}:{commit_hash}` to `{module_dir.as_posix()}`"
        )
        repo = git.Repo.clone_from(repository_url, module_dir.as_posix())
        repo.git.checkout(commit_hash)
    else:
        repo = git.Repo(module_dir.as_posix())

    if repo.head.commit.hexsha[:7] != commit_hash:
        logging.error(
            f"ERROR: Failed while cloning module `{repository_url}:{commit_hash}`"
        )
        logging.error(f"{commit_hash} does not match {repo.head.commit.hexsha[:7]}`")
        raise RuntimeError(
            f"ERROR: {commit_hash} does not match {repo.head.commit.hexsha[:7]}"
        )

    return module_dir


def dump_parameters_to_file(output_dir: Path, parameters: List[str]) -> None:
    os.makedirs(output_dir, exist_ok=True)

    if parameters is not None:
        params_file = os.path.join(output_dir, "parameters.txt")
        with open(params_file, "w") as params_file:
            params_file.write(f"{parameters}")

        param_dict_file = os.path.join(output_dir, "..", "parameters_dict.txt")
        with open(param_dict_file, "a") as param_dict_file:
            param_dict_file.write(f"{os.path.basename(output_dir)} {parameters}\n")


def post_execution(
    module_name: str,
    is_metric: bool,
    output_dir: Path,
    output_files: List[str],
    exit_code: int,
) -> None:
    if is_metric and len(output_files) > 0:
        for output_file_path in output_files:
            # File exist checks should be handled by snakemake already
            is_success = exit_code == 0  # and os.path.exists(output_file_path)
            output_file = os.path.basename(output_file_path)
            append_metric_mapping(module_name, output_dir, output_file, is_success)


def append_metric_mapping(
    metric_name: str, output_dir: Path, output_file: str, is_success: bool
) -> None:
    # Remove cwd from output_dir if it's a prefix
    cwd = Path(os.getcwd())
    if output_dir.is_relative_to(cwd):
        output_dir = output_dir.relative_to(cwd)

    # Split output_dir into root_output_dir and output_dir
    root_output_dir = output_dir.parts[0]
    output_dir = Path(*output_dir.parts[1:])

    metrics_mapping_file = Path(root_output_dir) / "metrics.mapping.tsv"

    metric_id = generate_metric_id(metric_name, output_dir, output_file)
    timestamp = datetime.now().isoformat()

    new_row = [metric_id, timestamp, metric_name, output_dir, output_file, is_success]

    file_exists = metrics_mapping_file.exists()
    with open(metrics_mapping_file, mode="a", newline="") as f:
        writer = csv.writer(f, delimiter="\t")

        if not file_exists:
            writer.writerow(
                [
                    "metric_id",
                    "timestamp",
                    "metric_name",
                    "output_path",
                    "output_file",
                    "is_success",
                ]
            )

        # Append the new row
        writer.writerow(new_row)


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

    keep_module_logs = params.get("keep_module_logs", False)

    output_dir = Path(str(os.path.commonpath(snakemake.output)))
    if len(snakemake.output) == 1:
        output_dir = Path(os.path.dirname(output_dir))

    # Create parameters file for outputs
    dump_parameters_to_file(output_dir, parameters)

    # Clone git repository
    repositories_dir = Path(".snakemake") / "repos"
    module_dir = clone_module(repositories_dir, repository_url, commit_hash)

    # Execute module code
    module_name = get_module_name_from_rule_name(snakemake.rule)

    exit_code = execution(
        module_dir,
        module_name=module_name,
        output_dir=output_dir,
        dataset=dataset,
        inputs_map=inputs_map,
        parameters=parameters,
        keep_module_logs=keep_module_logs,
    )

    post_execution(
        module_name=module_name,
        is_metric=True,
        output_dir=output_dir,
        output_files=snakemake.output,
        exit_code=exit_code,
    )


except NameError:
    raise RuntimeError(f"This script must be run from within a Snakemake workflow")
