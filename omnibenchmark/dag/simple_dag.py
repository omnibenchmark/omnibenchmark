"""Simple DAG implementation to replace NetworkX dependency."""

import heapq
from typing import Dict, List, Set, Tuple, Any, Iterator
from collections import defaultdict


class CyclicDependencyError(Exception):
    """Raised when a cycle is detected in the DAG."""

    pass


class SimpleDAG:
    """A simple directed acyclic graph implementation.

    Deterministic: nodes and edges keep insertion order (dicts used as ordered
    sets, never `set`, whose iteration order follows string hashes and so
    changes with PYTHONHASHSEED). Every iteration, and `topological_sort`'s tie
    order, is therefore the same on every run and every machine.
    """

    def __init__(self) -> None:
        """Initialize an empty DAG."""
        self.nodes: Dict[Any, None] = {}
        self._edges: Dict[Any, Dict[Any, None]] = defaultdict(dict)
        self.predecessors: Dict[Any, Dict[Any, None]] = defaultdict(dict)
        self.node_attrs: Dict[Any, Dict[str, Any]] = defaultdict(dict)

    def add_node(self, node: Any, **attrs: Any) -> None:
        """Add a node to the graph with optional attributes."""
        self.nodes.setdefault(node)
        self.node_attrs[node].update(attrs)

    def add_nodes_from(
        self, nodes_with_attrs: List[Tuple[Any, Dict[str, Any]]]
    ) -> None:
        """Add multiple nodes with attributes."""
        for node, attrs in nodes_with_attrs:
            self.add_node(node, **attrs)

    def add_edge(self, from_node: Any, to_node: Any) -> None:
        """Add an edge from from_node to to_node."""
        # Ensure both nodes exist
        self.nodes.setdefault(from_node)
        self.nodes.setdefault(to_node)

        # Add the edge
        self._edges[from_node].setdefault(to_node)
        self.predecessors[to_node].setdefault(from_node)

    def in_degree(self) -> Iterator[Tuple[Any, int]]:
        """Return an iterator of (node, in_degree) pairs."""
        for node in self.nodes:
            yield (node, len(self.predecessors[node]))

    def out_degree(self) -> Iterator[Tuple[Any, int]]:
        """Return an iterator of (node, out_degree) pairs."""
        for node in self.nodes:
            yield (node, len(self._edges[node]))

    def get_node_attributes(self, name: str, default: Any = None) -> Dict[Any, Any]:
        """Get a specific attribute for all nodes."""
        result: Dict[Any, Any] = {}
        for node in self.nodes:
            result[node] = self.node_attrs[node].get(name, default)
        return result

    def topological_sort(self) -> List[Any]:
        """
        Return a list of nodes in topological order.

        Stable: among the nodes whose predecessors are all placed, the one added
        first comes next. So the order is the same on every run, and a graph
        whose insertion order is already topological sorts to exactly that order.

        Raises:
            CyclicDependencyError: If the graph contains a cycle.
        """
        index = {node: i for i, node in enumerate(self.nodes)}
        in_degree = {node: len(self.predecessors[node]) for node in self.nodes}

        # Min-heap of insertion indices of the nodes ready to be placed.
        ready = [index[node] for node in self.nodes if in_degree[node] == 0]
        heapq.heapify(ready)
        order = list(self.nodes)
        result: List[Any] = []

        while ready:
            node = order[heapq.heappop(ready)]
            result.append(node)

            for neighbor in self._edges[node]:
                in_degree[neighbor] -= 1
                if in_degree[neighbor] == 0:
                    heapq.heappush(ready, index[neighbor])

        # If we haven't processed all nodes, there's a cycle
        if len(result) != len(self.nodes):
            raise CyclicDependencyError(
                "The graph contains a cycle and cannot be topologically sorted"
            )

        return result

    def all_simple_paths(self, source: Any, target: Any) -> List[List[Any]]:
        """Find all simple paths from source to target."""
        if source not in self.nodes or target not in self.nodes:
            return []

        paths: List[List[Any]] = []

        def dfs(current: Any, target: Any, path: List[Any], visited: Set[Any]) -> None:
            if current == target:
                paths.append(path[:])
                return

            for neighbor in self._edges[current]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    path.append(neighbor)
                    dfs(neighbor, target, path, visited)
                    path.pop()
                    visited.remove(neighbor)

        visited: Set[Any] = {source}
        dfs(source, target, [source], visited)
        return paths

    @property
    def edges(self) -> Iterator[Tuple[Any, Any]]:
        """Return an iterator of edge tuples (source, target)."""
        for source, targets in self._edges.items():
            for target in targets:
                yield (source, target)


# Compatibility exports
DiGraph = SimpleDAG
NetworkXUnfeasible = CyclicDependencyError


def topological_sort(graph: SimpleDAG) -> List[Any]:
    """Compatibility function for nx.topological_sort."""
    return graph.topological_sort()


def all_simple_paths(graph: SimpleDAG, source: Any, target: Any) -> List[List[Any]]:
    """Compatibility function for nx.all_simple_paths."""
    return graph.all_simple_paths(source, target)


def get_node_attributes(
    graph: SimpleDAG, name: str, default: Any = None
) -> Dict[Any, Any]:
    """Compatibility function for nx.get_node_attributes."""
    return graph.get_node_attributes(name, default)
