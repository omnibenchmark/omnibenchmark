"""
Backend module for omnibenchmark execution.

This module provides the compilation pipeline for benchmark execution:
- Module resolution (cloning, entrypoint detection)
- Workflow generation (Snakefile creation)
- Execution backends (Snakemake, etc.)
"""

from .resolver import ModuleResolver, ResolvedModule
from .cache import ModuleCache

__all__ = ["ModuleResolver", "ResolvedModule", "ModuleCache"]
