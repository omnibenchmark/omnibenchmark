# DAG Module

A lightweight directed acyclic graph implementation to replace NetworkX dependency.

## Features

- Basic DAG operations (add nodes/edges, topological sort, cycle detection)
- Path finding (all simple paths between nodes)
- Node attributes support
- NetworkX-compatible API for drop-in replacement

## Usage

```python
from omnibenchmark.dag import SimpleDAG

# Create a DAG
dag = SimpleDAG()

# Add nodes with attributes
dag.add_node("A", color="red")
dag.add_node("B", color="blue")

# Add edges
dag.add_edge("A", "B")
dag.add_edge("B", "C")

# Topological sort
order = dag.topological_sort()  # ['A', 'B', 'C']

# Find paths
paths = dag.all_simple_paths("A", "C")  # [['A', 'B', 'C']]

# Check degrees
for node, degree in dag.in_degree():
    print(f"{node}: {degree}")
```

## NetworkX Compatibility

The module provides drop-in replacements for common NetworkX operations:

- `DiGraph` → `SimpleDAG`
- `nx.topological_sort()` → `topological_sort()`
- `nx.all_simple_paths()` → `all_simple_paths()`
- `nx.NetworkXUnfeasible` → `CyclicDependencyError`
