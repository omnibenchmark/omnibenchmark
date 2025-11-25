"""Tests for SimpleDAG implementation."""

import pytest
from omnibenchmark.dag import (
    SimpleDAG,
    DiGraph,
    CyclicDependencyError,
    topological_sort,
    all_simple_paths,
    get_node_attributes,
)


@pytest.mark.short
class TestSimpleDAG:
    """Test the SimpleDAG implementation."""

    def test_initialization(self):
        """Test DAG initialization."""
        dag = SimpleDAG()
        assert len(dag.nodes) == 0
        assert len(list(dag.edges)) == 0
        assert len(dag._edges) == 0
        assert len(dag.predecessors) == 0

    def test_add_node(self):
        """Test adding nodes."""
        dag = SimpleDAG()
        dag.add_node("A")
        dag.add_node("B", color="red", weight=5)

        assert "A" in dag.nodes
        assert "B" in dag.nodes
        assert dag.node_attrs["B"]["color"] == "red"
        assert dag.node_attrs["B"]["weight"] == 5

    def test_add_nodes_from(self):
        """Test adding multiple nodes with attributes."""
        dag = SimpleDAG()
        nodes = [
            ("A", {"color": "red"}),
            ("B", {"color": "blue"}),
            ("C", {"color": "green"}),
        ]
        dag.add_nodes_from(nodes)

        assert len(dag.nodes) == 3
        assert dag.node_attrs["A"]["color"] == "red"
        assert dag.node_attrs["B"]["color"] == "blue"
        assert dag.node_attrs["C"]["color"] == "green"

    def test_add_edge(self):
        """Test adding edges."""
        dag = SimpleDAG()
        dag.add_edge("A", "B")
        dag.add_edge("B", "C")

        assert "B" in dag._edges["A"]
        assert "C" in dag._edges["B"]
        assert "A" in dag.predecessors["B"]
        assert "B" in dag.predecessors["C"]

    def test_in_degree(self):
        """Test in_degree calculation."""
        dag = SimpleDAG()
        dag.add_edge("A", "B")
        dag.add_edge("C", "B")
        dag.add_edge("B", "D")

        in_degrees = dict(dag.in_degree())
        assert in_degrees["A"] == 0
        assert in_degrees["B"] == 2
        assert in_degrees["C"] == 0
        assert in_degrees["D"] == 1

    def test_out_degree(self):
        """Test out_degree calculation."""
        dag = SimpleDAG()
        dag.add_edge("A", "B")
        dag.add_edge("A", "C")
        dag.add_edge("B", "D")

        out_degrees = dict(dag.out_degree())
        assert out_degrees["A"] == 2
        assert out_degrees["B"] == 1
        assert out_degrees["C"] == 0
        assert out_degrees["D"] == 0

    def test_get_node_attributes(self):
        """Test getting node attributes."""
        dag = SimpleDAG()
        dag.add_node("A", stage="input")
        dag.add_node("B", stage="process")
        dag.add_node("C", stage="output")
        dag.add_node("D")  # No stage attribute

        stages = dag.get_node_attributes("stage", default="none")
        assert stages["A"] == "input"
        assert stages["B"] == "process"
        assert stages["C"] == "output"
        assert stages["D"] == "none"

    def test_topological_sort_simple(self):
        """Test topological sort on a simple DAG."""
        dag = SimpleDAG()
        dag.add_edge("A", "B")
        dag.add_edge("B", "C")
        dag.add_edge("A", "C")

        order = dag.topological_sort()
        assert order.index("A") < order.index("B")
        assert order.index("B") < order.index("C")
        assert order.index("A") < order.index("C")

    def test_topological_sort_complex(self):
        """Test topological sort on a more complex DAG."""
        dag = SimpleDAG()
        # Create a diamond pattern
        dag.add_edge("A", "B")
        dag.add_edge("A", "C")
        dag.add_edge("B", "D")
        dag.add_edge("C", "D")

        order = dag.topological_sort()
        assert order.index("A") < order.index("B")
        assert order.index("A") < order.index("C")
        assert order.index("B") < order.index("D")
        assert order.index("C") < order.index("D")

    def test_topological_sort_with_cycle(self):
        """Test that topological sort detects cycles."""
        dag = SimpleDAG()
        dag.add_edge("A", "B")
        dag.add_edge("B", "C")
        dag.add_edge("C", "A")  # Creates a cycle

        with pytest.raises(CyclicDependencyError):
            dag.topological_sort()

    def test_all_simple_paths_basic(self):
        """Test finding all simple paths."""
        dag = SimpleDAG()
        dag.add_edge("A", "B")
        dag.add_edge("B", "C")
        dag.add_edge("A", "C")

        paths = dag.all_simple_paths("A", "C")
        assert len(paths) == 2
        assert ["A", "C"] in paths
        assert ["A", "B", "C"] in paths

    def test_all_simple_paths_diamond(self):
        """Test finding paths in diamond pattern."""
        dag = SimpleDAG()
        dag.add_edge("A", "B")
        dag.add_edge("A", "C")
        dag.add_edge("B", "D")
        dag.add_edge("C", "D")

        paths = dag.all_simple_paths("A", "D")
        assert len(paths) == 2
        assert ["A", "B", "D"] in paths
        assert ["A", "C", "D"] in paths

    def test_all_simple_paths_no_path(self):
        """Test when there's no path between nodes."""
        dag = SimpleDAG()
        dag.add_edge("A", "B")
        dag.add_edge("C", "D")

        paths = dag.all_simple_paths("A", "D")
        assert len(paths) == 0

    def test_all_simple_paths_nonexistent_nodes(self):
        """Test paths with nonexistent nodes."""
        dag = SimpleDAG()
        dag.add_edge("A", "B")

        paths = dag.all_simple_paths("A", "Z")
        assert len(paths) == 0

        paths = dag.all_simple_paths("Z", "B")
        assert len(paths) == 0

    def test_compatibility_functions(self):
        """Test the compatibility functions work as expected."""
        dag = SimpleDAG()
        dag.add_edge("A", "B")
        dag.add_edge("B", "C")
        dag.add_node("D", color="red")

        # Test topological_sort function
        order = topological_sort(dag)
        assert order.index("A") < order.index("B") < order.index("C")

        # Test all_simple_paths function
        paths = all_simple_paths(dag, "A", "C")
        assert len(paths) == 1
        assert paths[0] == ["A", "B", "C"]

        # Test get_node_attributes function
        colors = get_node_attributes(dag, "color", default="none")
        assert colors["D"] == "red"
        assert colors["A"] == "none"

    def test_diraph_alias(self):
        """Test that DiGraph is an alias for SimpleDAG."""
        dag = DiGraph()
        dag.add_edge("A", "B")
        assert "B" in dag._edges["A"]
        assert isinstance(dag, SimpleDAG)

    def test_isolated_nodes(self):
        """Test handling of isolated nodes."""
        dag = SimpleDAG()
        dag.add_node("A")
        dag.add_node("B")
        dag.add_edge("C", "D")

        # Isolated nodes should appear in topological sort
        order = dag.topological_sort()
        assert "A" in order
        assert "B" in order
        assert "C" in order
        assert "D" in order
        assert order.index("C") < order.index("D")

    def test_self_loop_detection(self):
        """Test that self-loops are detected as cycles."""
        dag = SimpleDAG()
        dag.add_edge("A", "A")  # Self-loop

        with pytest.raises(CyclicDependencyError):
            dag.topological_sort()

    def test_multiple_paths_complex(self):
        """Test finding paths in a more complex graph."""
        dag = SimpleDAG()
        # Create a graph with multiple paths from A to E
        dag.add_edge("A", "B")
        dag.add_edge("A", "C")
        dag.add_edge("B", "D")
        dag.add_edge("C", "D")
        dag.add_edge("B", "E")
        dag.add_edge("D", "E")

        paths = dag.all_simple_paths("A", "E")
        assert len(paths) == 3
        # Check that all expected paths are found
        expected_paths = [["A", "B", "E"], ["A", "B", "D", "E"], ["A", "C", "D", "E"]]
        for path in expected_paths:
            assert path in paths

    def test_edges_property(self):
        """Test that edges property returns correct edge tuples."""
        dag = SimpleDAG()
        dag.add_edge("A", "B")
        dag.add_edge("A", "C")
        dag.add_edge("B", "D")

        edges = list(dag.edges)
        assert len(edges) == 3
        assert ("A", "B") in edges
        assert ("A", "C") in edges
        assert ("B", "D") in edges
