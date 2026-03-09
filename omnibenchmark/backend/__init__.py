"""
Backend module for omnibenchmark execution.

Provides the compilation pipeline:
- Module resolution (cloning, entrypoint detection)
- Snakefile generation from resolved nodes
"""

from .resolver import ModuleResolver

__all__ = ["ModuleResolver"]
