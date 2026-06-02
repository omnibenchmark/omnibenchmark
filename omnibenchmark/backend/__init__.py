"""
Execution backend — compile a benchmark into a runnable workflow.

Represents: the compilation pipeline that turns resolved nodes into an
executable Snakemake workflow:
- Module resolution (cloning, entrypoint detection)
- Snakefile generation from resolved nodes
Layer: service
Depends on: core, model, git.
"""

from .resolver import ModuleResolver

__all__ = ["ModuleResolver"]
