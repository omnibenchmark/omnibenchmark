#! /usr/bin/env python
# WARNING: Custom dependencies might not be available here, since this is run inside a specified environment.

# CRITICAL: Modify sys.path FIRST, before any other imports
# In Snakemake 9+, when running in containers, the omnibenchmark package root
# may be mounted but not added to sys.path. We need to add it before importing.
import sys
from pathlib import Path

_script_dir = Path(__file__).parent
_omnibenchmark_root = _script_dir.parent.parent.parent.parent
if str(_omnibenchmark_root) not in sys.path:
    sys.path.insert(0, str(_omnibenchmark_root))

# Now safe to import everything else
import logging
import os
import traceback
from typing import Any, Dict, Optional

# ruff: noqa: E402
# Imports below must come after sys.path modification above
from snakemake.script import Snakemake

from omnibenchmark.benchmark import constants
from omnibenchmark.benchmark.params import Params
from omnibenchmark.benchmark.symlinks import SymlinkManager
from omnibenchmark.git.clone import clone_module
from omnibenchmark.workflow.snakemake.scripts.execution import execution
from omnibenchmark.workflow.snakemake.scripts.inputs_map_builder import build_inputs_map

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
# Use attribute access for inputs_map to allow Snakemake to evaluate lambda functions
inputs_map_template: Optional[Dict[str, Any]] = getattr(
    snakemake.params, "inputs_map", None
)
dataset: str = params.get("dataset", getattr(snakemake.wildcards, "dataset", "unknown"))

# Build inputs_map using actual file paths from snakemake.input
# This is critical for remote storage where files are downloaded to .snakemake/storage/
inputs_map: Dict[str, Any] = {}
if inputs_map_template:
    inputs_map = build_inputs_map(inputs_map_template, list(snakemake.input))

# For now we're handling timeout in seconds.
# When implementing cluster resource handling, we needt to convert this to minutes (e.g. slurm takes it in min)
timeout: Optional[int] = params.get(constants.LOCAL_TIMEOUT_VAR, None)

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
if repository_url is None or commit_hash is None:
    raise RuntimeError("repository_url and commit_hash must be provided")
module_dir = clone_module(repository_url, commit_hash)

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
