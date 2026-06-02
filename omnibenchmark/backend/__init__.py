"""
Execution backend — compile a benchmark into a runnable workflow.

Represents: the compilation pipeline that turns resolved nodes into an
executable Snakemake workflow.
Layer: service
Depends on: core, model, git.

Pipeline: module resolution (cloning, entrypoint detection), then Snakefile
generation from the resolved nodes.
"""

from .resolver import ModuleResolver

__all__ = ["ModuleResolver"]
