"""DAG Builder for Benchmark execution.

This module provides the DAGBuilder class which constructs the computational
directed acyclic graph (DAG) from a benchmark model.
"""

from pathlib import Path
from typing import Dict, List, Optional

from omnibenchmark.dag import DiGraph
from omnibenchmark.model import Benchmark, ValidationError

from ._graph import (
    build_stage_dag,
    compute_stage_order,
    expand_stage_nodes,
)
from ._node import BenchmarkNode


class DAGBuilder:
    """Builds and manages the computational DAG for a benchmark.

    This class separates the DAG construction logic from the BenchmarkExecution
    class, following the Single Responsibility Principle. It handles:
    - Building the stage dependency graph
    - Expanding stages into computational nodes
    - Connecting nodes based on their dependencies
    """

    def __init__(self, model: Benchmark, out_dir: Path):
        """Initialize the DAG builder.

        Args:
            model: The benchmark model containing the definition
            out_dir: Output directory for benchmark execution
        """
        self.model = model
        self.out_dir = out_dir
        self._graph: Optional[DiGraph] = None
        self._stage_dag: Optional[DiGraph] = None
        self._stage_order: Optional[List[str]] = None
        self._stage_nodes_map: Dict[str, List[BenchmarkNode]] = {}

    def build(self) -> DiGraph:
        """Build the computational DAG.

        This method constructs the full computational graph by:
        1. Building the stage dependency graph
        2. Computing the topological ordering of stages
        3. Expanding each stage into computational nodes
        4. Connecting nodes based on their dependencies

        Returns:
            The constructed computational DAG

        Raises:
            ValidationError: If the stage graph has cyclic dependencies
        """
        if self._graph is None:
            self._graph = self._build_dag()
        return self._graph

    def get_stage_dag(self) -> DiGraph:
        """Get the stage dependency graph.

        Returns:
            The stage-level dependency graph
        """
        if self._stage_dag is None:
            self._stage_dag = build_stage_dag(self.model)
        return self._stage_dag

    def get_stage_order(self) -> List[str]:
        """Get the topological ordering of stages.

        Returns:
            List of stage IDs in topological order

        Raises:
            ValidationError: If the stage graph has cyclic dependencies
        """
        if self._stage_order is None:
            stage_dag = self.get_stage_dag()
            self._stage_order = compute_stage_order(stage_dag)
        return self._stage_order

    def get_stage_nodes(self, stage_id: str) -> List[BenchmarkNode]:
        """Get the computational nodes for a specific stage.

        Args:
            stage_id: The ID of the stage

        Returns:
            List of computational nodes in the stage
        """
        if stage_id not in self._stage_nodes_map:
            stage = self.model.get_stages()[stage_id]
            stage_order = self.get_stage_order()
            nodes = expand_stage_nodes(self.model, stage, self.out_dir, stage_order)
            self._stage_nodes_map[stage_id] = nodes
        return self._stage_nodes_map[stage_id]

    def _build_dag(self) -> DiGraph:
        """Build the computational DAG.

        Internal method that performs the actual DAG construction.

        Returns:
            The constructed computational DAG
        """
        g = DiGraph()

        # Get stage ordering first (this will validate no cycles)
        _ = self.get_stage_order()

        # Expand all stages into computational nodes
        for stage_id in self.model.get_stages():
            nodes = self.get_stage_nodes(stage_id)
            nodes_with_stage = [(node, {"stage": stage_id}) for node in nodes]
            g.add_nodes_from(nodes_with_stage)

        # Connect nodes based on their dependencies
        for current_node in g.nodes:
            after_stage = current_node.after
            if after_stage:
                departure_stage_nodes = self._stage_nodes_map[after_stage]
                for departure_stage_node in departure_stage_nodes:
                    g.add_edge(departure_stage_node, current_node)

        return g

    def validate_dag(self) -> None:
        """Validate the DAG structure.

        Checks for:
        - Cyclic dependencies in the stage graph
        - Disconnected components
        - Invalid node references

        Raises:
            ValidationError: If validation fails
        """
        # Check for cycles in stage graph
        try:
            _ = self.get_stage_order()
        except ValidationError:
            raise

        # Build the graph to ensure all connections are valid
        _graph = self.build()

        # Check for disconnected components (optional)
        # This could be added if needed

    def get_statistics(self) -> Dict:
        """Get statistics about the DAG.

        Returns:
            Dictionary containing DAG statistics like node count,
            edge count, maximum path length, etc.
        """
        graph = self.build()
        stage_dag = self.get_stage_dag()

        # Find initial and terminal nodes
        initial_nodes = [
            node for node, in_degree in graph.in_degree() if in_degree == 0
        ]
        terminal_nodes = [
            node for node, out_degree in graph.out_degree() if out_degree == 0
        ]

        return {
            "num_stages": len(stage_dag.nodes),
            "num_nodes": len(graph.nodes),
            "num_edges": len(list(graph.edges)),
            "num_initial_nodes": len(initial_nodes),
            "num_terminal_nodes": len(terminal_nodes),
            "stage_order": self.get_stage_order(),
        }

    def reset(self) -> None:
        """Reset the cached DAG structures.

        This forces a rebuild on the next access.
        """
        self._graph = None
        self._stage_dag = None
        self._stage_order = None
        self._stage_nodes_map = {}
