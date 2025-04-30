#! /usr/bin/env python
# WARNING: Custom dependencies might not be available here, since this is run inside a specified environment.

import os
from pathlib import Path

from snakemake.script import Snakemake

from omnibenchmark.benchmark import constants
from omnibenchmark.io.code import clone_module
from omnibenchmark.workflow.snakemake.scripts.execution import execution
from omnibenchmark.benchmark.symlinks import SymlinkManager
from omnibenchmark.workflow.snakemake.scripts.utils import (
    get_module_name_from_rule_name,
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

    manager = SymlinkManager(output_dir.parent)
    manager.store(parameters)

    # Clone git repository
    repositories_dir = Path(".snakemake") / "repos"
    module_dir = clone_module(repositories_dir, repository_url, commit_hash)

    # Execute module code
    module_name = get_module_name_from_rule_name(snakemake.rule)

    resources = dict(snakemake.resources) if hasattr(snakemake, "resources") else {}
    # For now we're handling timeout in seconds as runtime.
    # When implementing cluster resource handling, we needt to convert this to minutes (e.g. slurm takes it in min)
    timeout = resources.get("runtime", constants.DEFAULT_TIMEOUT_SECONDS)

    # TODO(ben): we don't do anything with this exit_code?
    exit_code = execution(
        module_dir,
        module_name=module_name,
        output_dir=output_dir,
        dataset=dataset,
        inputs_map=inputs_map,
        parameters=parameters,
        keep_module_logs=keep_module_logs,
        timeout=timeout,
    )


except NameError:
    raise RuntimeError("This script must be run from within a Snakemake workflow")
