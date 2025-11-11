"""Tests for the DAGBuilder class."""

import pytest
from unittest.mock import Mock

from omnibenchmark.benchmark._dag_builder import DAGBuilder
from omnibenchmark.benchmark.params import Params
from omnibenchmark.dag import DiGraph
from omnibenchmark.model import Benchmark, Stage, ValidationError


@pytest.fixture
def simple_benchmark():
    """Create a simple benchmark model for testing."""
    # Create a mock benchmark with minimal structure
    benchmark = Mock(spec=Benchmark)

    # Create stages
    stage1 = Mock(spec=Stage)
    stage1.id = "data"
    stage1.inputs = None
    stage1.modules = [Mock(id="dataset1")]

    stage2 = Mock(spec=Stage)
    stage2.id = "process"
    stage2.inputs = [Mock(entries=["data_output"])]
    stage2.modules = [Mock(id="method1")]

    stages = {"data": stage1, "process": stage2}

    # Setup benchmark methods
    benchmark.get_stages.return_value = stages
    benchmark.is_initial.side_effect = lambda s: s.id == "data"
    benchmark.get_stage_outputs.return_value = {"output": "data_output"}
    benchmark.get_stage_implicit_inputs.side_effect = (
        lambda s: [] if s.id == "data" else [["data_output"]]
    )
    benchmark.get_modules_by_stage.side_effect = lambda s: {m.id: m for m in s.modules}
    benchmark.get_module_parameters.return_value = [Params({"param": "value"})]
    benchmark.get_output_stage.side_effect = (
        lambda x: stage1 if x == "data_output" else None
    )
    benchmark.get_stage_by_output.side_effect = (
        lambda x: stage1 if x == "data_output" else None
    )
    benchmark.get_explicit_inputs.return_value = {"input": "data_output"}

    return benchmark


@pytest.fixture
def out_dir(tmp_path):
    """Create a temporary output directory."""
    return tmp_path / "out"


class TestDAGBuilder:
    """Test suite for DAGBuilder class."""

    @pytest.mark.short
    def test_initialization(self, simple_benchmark, out_dir):
        """Test DAGBuilder initialization."""
        builder = DAGBuilder(simple_benchmark, out_dir)

        assert builder.model == simple_benchmark
        assert builder.out_dir == out_dir
        assert builder._graph is None
        assert builder._stage_dag is None
        assert builder._stage_order is None
        assert builder._stage_nodes_map == {}

    @pytest.mark.short
    def test_build_creates_graph(self, simple_benchmark, out_dir):
        """Test that build() creates a DiGraph."""
        builder = DAGBuilder(simple_benchmark, out_dir)
        graph = builder.build()

        assert isinstance(graph, DiGraph)
        assert builder._graph is not None
        assert builder._graph == graph

    @pytest.mark.short
    def test_build_caches_result(self, simple_benchmark, out_dir):
        """Test that build() caches and returns the same graph."""
        builder = DAGBuilder(simple_benchmark, out_dir)

        graph1 = builder.build()
        graph2 = builder.build()

        assert graph1 is graph2

    @pytest.mark.short
    def test_get_stage_dag(self, simple_benchmark, out_dir):
        """Test stage DAG construction."""
        builder = DAGBuilder(simple_benchmark, out_dir)
        stage_dag = builder.get_stage_dag()

        assert isinstance(stage_dag, DiGraph)
        assert "data" in stage_dag.nodes
        assert "process" in stage_dag.nodes

    @pytest.mark.short
    def test_get_stage_order(self, simple_benchmark, out_dir):
        """Test topological ordering of stages."""
        builder = DAGBuilder(simple_benchmark, out_dir)
        order = builder.get_stage_order()

        assert isinstance(order, list)
        assert "data" in order
        assert "process" in order
        # Initial stage should come before process stage
        assert order.index("data") < order.index("process")

    @pytest.mark.short
    def test_get_stage_nodes(self, simple_benchmark, out_dir):
        """Test getting nodes for a specific stage."""
        builder = DAGBuilder(simple_benchmark, out_dir)

        # Build the graph first to populate stage nodes
        builder.build()

        data_nodes = builder.get_stage_nodes("data")
        assert isinstance(data_nodes, list)
        assert len(data_nodes) > 0

        # Check caching
        data_nodes2 = builder.get_stage_nodes("data")
        assert data_nodes == data_nodes2

    @pytest.mark.short
    def test_reset(self, simple_benchmark, out_dir):
        """Test resetting the builder state."""
        builder = DAGBuilder(simple_benchmark, out_dir)

        # Build graph and populate caches
        _graph1 = builder.build()
        _ = builder.get_stage_order()

        # Reset
        builder.reset()

        assert builder._graph is None
        assert builder._stage_dag is None
        assert builder._stage_order is None
        assert builder._stage_nodes_map == {}

        # Build again should work
        graph2 = builder.build()
        assert isinstance(graph2, DiGraph)

    @pytest.mark.short
    def test_get_statistics(self, simple_benchmark, out_dir):
        """Test getting DAG statistics."""
        builder = DAGBuilder(simple_benchmark, out_dir)
        stats = builder.get_statistics()

        assert isinstance(stats, dict)
        assert "num_stages" in stats
        assert "num_nodes" in stats
        assert "num_edges" in stats
        assert "num_initial_nodes" in stats
        assert "num_terminal_nodes" in stats
        assert "stage_order" in stats

        assert stats["num_stages"] == 2  # data and process
        assert stats["num_nodes"] >= 2  # at least one node per stage
        assert isinstance(stats["stage_order"], list)

    @pytest.mark.short
    def test_validate_dag_with_valid_graph(self, simple_benchmark, out_dir):
        """Test DAG validation with a valid graph."""
        builder = DAGBuilder(simple_benchmark, out_dir)

        # Should not raise any exceptions
        builder.validate_dag()

    @pytest.mark.short
    def test_validate_dag_with_cycle(self, out_dir):
        """Test DAG validation with cyclic dependencies."""
        # Create a benchmark with cyclic stage dependencies
        benchmark = Mock(spec=Benchmark)

        stage1 = Mock(spec=Stage)
        stage1.id = "stage1"
        stage1.inputs = [Mock(entries=["output2"])]

        stage2 = Mock(spec=Stage)
        stage2.id = "stage2"
        stage2.inputs = [Mock(entries=["output1"])]

        benchmark.get_stages.return_value = {"stage1": stage1, "stage2": stage2}
        benchmark.get_output_stage.side_effect = (
            lambda x: stage1 if x == "output1" else stage2
        )

        builder = DAGBuilder(benchmark, out_dir)

        with pytest.raises(ValidationError, match="cyclic dependencies"):
            builder.validate_dag()

    @pytest.mark.short
    def test_empty_benchmark(self, out_dir):
        """Test building DAG with empty benchmark."""
        benchmark = Mock(spec=Benchmark)
        benchmark.get_stages.return_value = {}

        builder = DAGBuilder(benchmark, out_dir)
        graph = builder.build()

        assert isinstance(graph, DiGraph)
        assert len(graph.nodes) == 0
        assert len(list(graph.edges)) == 0

    @pytest.mark.short
    def test_single_stage_benchmark(self, out_dir):
        """Test building DAG with single stage."""
        benchmark = Mock(spec=Benchmark)

        stage = Mock(spec=Stage)
        stage.id = "single"
        stage.inputs = None
        stage.modules = [Mock(id="module1")]

        benchmark.get_stages.return_value = {"single": stage}
        benchmark.is_initial.return_value = True
        benchmark.get_stage_outputs.return_value = {"output": "single_output"}
        benchmark.get_stage_implicit_inputs.return_value = []
        benchmark.get_modules_by_stage.return_value = {"module1": stage.modules[0]}
        benchmark.get_module_parameters.return_value = None

        builder = DAGBuilder(benchmark, out_dir)
        graph = builder.build()

        assert isinstance(graph, DiGraph)
        assert len(graph.nodes) >= 1

        stats = builder.get_statistics()
        assert stats["num_stages"] == 1
        assert stats["num_initial_nodes"] >= 1
        assert stats["num_terminal_nodes"] >= 1
