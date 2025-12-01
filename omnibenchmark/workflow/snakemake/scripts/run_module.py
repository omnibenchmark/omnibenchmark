#! /usr/bin/env python
# WARNING: Custom dependencies might not be available here, since this is run inside a specified environment.

import logging
import os
import sys
import traceback
from typing import Any, Dict, Optional

from pathlib import Path

from snakemake.script import Snakemake

from omnibenchmark.benchmark import constants
from omnibenchmark.benchmark.params import Params
from omnibenchmark.io.code import clone_module
from omnibenchmark.workflow.snakemake.scripts.execution import execution
from omnibenchmark.benchmark.symlinks import SymlinkManager

logger = logging.getLogger("SNAKEMAKE_RUNNER")

try:
    snakemake: Snakemake = snakemake  # type: ignore[name-defined]
except NameError:
    raise RuntimeError("This script must be run from within a Snakemake workflow")

params: Dict[str, Any] = dict(snakemake.params)

def extract_module_name(rule_name: str) -> str:
    parts = rule_name.split("_")
    if len(parts) == 0:
        raise ValueError("Invalid rule name")
    return parts[1]



repository_url: Optional[str] = params.get("repository_url")
commit_hash: Optional[str] = params.get("commit_hash")
parameters: Optional[Any] = params.get("parameters")
inputs_map: Optional[Dict[str, Any]] = params.get("inputs_map")
dataset: str = params.get("dataset", getattr(snakemake.wildcards, "dataset", "unknown"))

# For now we're handling timeout in seconds.
# When implementing cluster resource handling, we needt to convert this to minutes (e.g. slurm takes it in min)
timeout: int = params.get(
    constants.LOCAL_TIMEOUT_VAR, constants.DEFAULT_TIMEOUT_SECONDS
)

keep_module_logs: bool = params.get("keep_module_logs", False)
keep_going: bool = snakemake.config.get("keep_going", False)

output_dir = Path(str(os.path.commonpath(snakemake.output)))
if len(snakemake.output) == 1:
    output_dir = Path(os.path.dirname(output_dir))

manager = SymlinkManager(output_dir.parent)
symlink_info = manager.store(parameters)

# Use symlink path instead of original output_dir
if symlink_info:
    output_dir = symlink_info["symlink_path"]

# Clone git repository
repositories_dir = Path(".snakemake") / "repos"
if repository_url is None or commit_hash is None:
    raise RuntimeError("repository_url and commit_hash must be provided")
module_dir = clone_module(repositories_dir, repository_url, commit_hash)

# Execute module code
module_name = extract_module_name(snakemake.rule)

try:
    # Handle None parameters and inputs_map by providing defaults
    if inputs_map is None:
        inputs_map = {}
    if parameters is None:
        parameters = Params()

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
