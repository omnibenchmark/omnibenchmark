#! /usr/bin/env python
# WARNING: Custom dependencies might not be available here, since this is run inside a specified environment.

import logging
import os
from pathlib import Path

from snakemake.script import Snakemake

from omni.io.code import clone_module
from omni.workflow.snakemake.scripts.execution import execution
from omni.workflow.snakemake.scripts.utils import (
    get_module_name_from_rule_name,
    dump_parameters_to_file,
)


try:
    snakemake: Snakemake = snakemake
    params = dict(snakemake.params)

    repository_url = params["repository_url"]
    commit_hash = params["commit_hash"]
    parameters = params.get("parameters")
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


except NameError:
    raise RuntimeError(f"This script must be run from within a Snakemake workflow")
