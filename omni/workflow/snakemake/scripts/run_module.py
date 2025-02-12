#! /usr/bin/env python
# WARNING: Custom dependencies might not be available here, since this is run inside a specified environment.

import os
from pathlib import Path

from snakemake.script import Snakemake

from omni.io.code import clone_module, dump_parameters_to_file
from omni.workflow.snakemake.scripts.execution import execution

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
        keep_module_logs=keep_module_logs,
    )

except NameError:
    raise RuntimeError(f"This script must be run from within a Snakemake workflow")
