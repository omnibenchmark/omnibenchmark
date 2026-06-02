"""
Core execution engine for omnibenchmark.

Represents: the engine that turns a validated benchmark definition into a run —
builds the DAG of stages/modules and drives it to completion.
Layer: engine (central domain)
Depends on: model, git, dag.

This is the central domain package: it turns a validated benchmark
definition (see :mod:`omnibenchmark.model`) into an executable workflow and
drives it to completion.

Key responsibilities:
- ``execution.py``  — :class:`BenchmarkExecution`, the top-level orchestrator
- ``_node.py`` / ``_graph.py`` / ``_dag_builder.py`` — build and traverse the
  benchmark DAG of stages and modules
- ``prefetch.py``   — populate the git cache with every repo a benchmark needs
- ``validate.py`` / ``validator.py`` — run-time validation of a benchmark
- ``cite.py`` / ``metadata.py`` — citation and provenance metadata
- ``dashboard.py`` / ``_dot.py`` / ``_mermaid.py`` — reporting and visualisation

Layering: ``core`` depends downward on :mod:`omnibenchmark.model`,
:mod:`omnibenchmark.git`, and the leaf utilities (``logging``, ``progress``,
``formatting``). Nothing in a lower layer imports ``core``.

(Formerly ``omnibenchmark.benchmark``; renamed to avoid overloading the word
"benchmark" across the tool, package, module, and model class.)
"""

from .execution import BenchmarkExecution
from .execution import BenchmarkExecution as Benchmark  # Compatibility alias
from ._node import BenchmarkNode
from .validator import Validator
from .prefetch import populate_git_cache


__all__ = [
    "BenchmarkExecution",
    "Benchmark",
    "BenchmarkNode",
    "Validator",
    "populate_git_cache",
]
