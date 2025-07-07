#! /usr/bin/env python
# WARNING: Custom dependencies might not be available here, since this is run inside a specified environment.

import logging
import os
import sys
import traceback

from pathlib import Path

from snakemake.script import Snakemake

from omnibenchmark.benchmark import constants
from omnibenchmark.io.code import clone_module
from omnibenchmark.workflow.snakemake.scripts.execution import execution
from omnibenchmark.benchmark.symlinks import SymlinkManager
from omnibenchmark.workflow.snakemake.scripts.utils import (
    get_module_name_from_rule_name,
)

logger = logging.getLogger("SNAKEMAKE_RUNNER")

try:
    snakemake: Snakemake = snakemake
except NameError:
    raise RuntimeError("This script must be run from within a Snakemake workflow")

params = dict(snakemake.params)

repository_url = params.get("repository_url")
commit_hash = params.get("commit_hash")
parameters = params.get("parameters")
inputs_map = params.get("inputs_map")
dataset = params.get("dataset", getattr(snakemake.wildcards, "dataset", "unknown"))

# For now we're handling timeout in seconds.
# When implementing cluster resource handling, we needt to convert this to minutes (e.g. slurm takes it in min)
timeout = params.get(constants.LOCAL_TIMEOUT_VAR, constants.DEFAULT_TIMEOUT_SECONDS)

keep_module_logs = params.get("keep_module_logs", False)
keep_going = snakemake.config.get("keep_going", False)

# For now we're handling timeout in seconds as runtime.
# When implementing cluster resource handling, we will need to convert this to minutes (e.g. slurm takes it in min)
timeout = params.get(constants.LOCAL_TIMEOUT_VAR, constants.DEFAULT_TIMEOUT_SECONDS)

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

try:
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
    if exit_code != 0:
        logger.warning(f"Module execution failed with exit code {exit_code}")
        if not keep_going:
            sys.exit(1)

except Exception as e:
    traceback.print_exc()
    logger.error(f"CRITICAL: {e}")
    if not keep_going:
        # explicit raise will leave a clearer traceback
        raise

sys.exit(0)
