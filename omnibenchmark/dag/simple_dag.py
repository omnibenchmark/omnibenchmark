"""Simple DAG implementation to replace NetworkX dependency."""

from typing import Dict, List, Set, Tuple, Any, Iterator
from collections import defaultdict, deque


class CyclicDependencyError(Exception):
    """Raised when a cycle is detected in the DAG."""

    pass


class SimpleDAG:
    """A simple directed acyclic graph implementation."""

    def __init__(self) -> None:
        """Initialize an empty DAG."""
        self.nodes: Set[Any] = set()
        self._edges: Dict[Any, Set[Any]] = defaultdict(set)
        self.predecessors: Dict[Any, Set[Any]] = defaultdict(set)
        self.node_attrs: Dict[Any, Dict[str, Any]] = defaultdict(dict)

    def add_node(self, node: Any, **attrs: Any) -> None:
        """Add a node to the graph with optional attributes."""
        self.nodes.add(node)
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
        self.nodes.add(from_node)
        self.nodes.add(to_node)

        # Add the edge
        self._edges[from_node].add(to_node)
        self.predecessors[to_node].add(from_node)

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

        Raises:
            CyclicDependencyError: If the graph contains a cycle.
        """
        # Count in-degrees
        in_degree: Dict[Any, int] = defaultdict(int)
        for node in self.nodes:
            in_degree[node] = len(self.predecessors[node])

        # Find all nodes with no incoming edges
        queue: deque[Any] = deque([node for node in self.nodes if in_degree[node] == 0])
        result: List[Any] = []

        while queue:
            node = queue.popleft()
            result.append(node)

            # For each neighbor, reduce its in-degree
            for neighbor in self._edges[node]:
                in_degree[neighbor] -= 1
                if in_degree[neighbor] == 0:
                    queue.append(neighbor)

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
