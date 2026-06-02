"""
Generic directed-acyclic-graph primitives.

Represents: backend-agnostic graph operations (topological sort, simple paths,
node attributes) — a small standalone toolkit, not benchmark-specific.
Layer: foundation (generic)
Depends on: nothing internal.
"""

from .simple_dag import (
    SimpleDAG,
    DiGraph,
    CyclicDependencyError,
    NetworkXUnfeasible,
    topological_sort,
    all_simple_paths,
    get_node_attributes,
)

__all__ = [
    "SimpleDAG",
    "DiGraph",
    "CyclicDependencyError",
    "NetworkXUnfeasible",
    "topological_sort",
    "all_simple_paths",
    "get_node_attributes",
]
