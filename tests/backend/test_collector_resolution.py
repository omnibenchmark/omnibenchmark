"""
Unit tests for metric collector resolution.

Tests the conversion of MetricCollector entities into ResolvedNode instances
through the collector_resolution module.
"""

import pytest
from pathlib import Path
from unittest.mock import Mock

from omnibenchmark.backend.collector_resolution import (
    resolve_metric_collectors,
    _expand_collector_parameters,
    _gather_collector_inputs,
    _resolve_collector_outputs,
)
from omnibenchmark.model.benchmark import MetricCollector, Repository, IOFile, Parameter
from omnibenchmark.model.resolved import (
    ResolvedNode,
    ResolvedModule,
    ResolvedEnvironment,
)
from omnibenchmark.model import SoftwareBackendEnum
from omnibenchmark.benchmark.params import Params


@pytest.fixture
def mock_resolver():
    """Create a mock ModuleResolver."""
    resolver = Mock()

    # Mock resolve method to return a ResolvedModule
    def mock_resolve(module, module_id, software_environment_id, dirty=False):
        return ResolvedModule(
            repository_url="https://github.com/test/collector.git",
            commit="abc123",
            module_dir=Path(".repos/abc123"),
            entrypoint=Path("collect.py"),
            software_environment_id=software_environment_id,
            resolved_environment=ResolvedEnvironment(
                backend_type=SoftwareBackendEnum.CONDA, reference="envs/collector.yaml"
            ),
            has_shebang=False,
            interpreter="python3",
        )

    resolver.resolve = Mock(side_effect=mock_resolve)
    return resolver


@pytest.fixture
def simple_collector():
    """Create a simple MetricCollector for testing."""
    return MetricCollector(
        id="plotting",
        name="Plotting collector",
        repository=Repository(
            url="https://github.com/test/collector.git", commit="abc123"
        ),
        software_environment="python",
        inputs=["metrics.scores"],
        outputs=[IOFile(id="report.html", path="{input}/{name}/report.html")],
    )


@pytest.fixture
def collector_with_params():
    """Create a collector with parameters."""
    return MetricCollector(
        id="analysis",
        name="Analysis collector",
        repository=Repository(
            url="https://github.com/test/analysis.git", commit="def456"
        ),
        software_environment="r",
        inputs=["results.data"],
        outputs=[IOFile(id="output.pdf", path="{name}/output.pdf")],
        parameters=[Parameter(id="threshold", values=[0.05, 0.01])],
    )


@pytest.fixture
def mock_benchmark():
    """Create a mock Benchmark with stages."""
    benchmark = Mock()

    # Mock stage with output
    mock_stage = Mock()
    mock_stage.id = "metrics"
    mock_stage.outputs = [Mock(id="scores")]

    def get_stage_by_output(output_id):
        if output_id == "metrics.scores" or output_id == "scores":
            return mock_stage
        elif output_id == "results.data" or output_id == "data":
            mock_stage2 = Mock()
            mock_stage2.id = "results"
            return mock_stage2
        return None

    benchmark.get_stage_by_output = Mock(side_effect=get_stage_by_output)
    benchmark.get_name = Mock(return_value="TestBenchmark")
    benchmark.get_version = Mock(return_value="1.0")
    benchmark.get_author = Mock(return_value="Test Author")

    return benchmark


@pytest.fixture
def sample_resolved_nodes():
    """Create sample resolved nodes that collectors will depend on."""
    module = ResolvedModule(
        repository_url="https://github.com/test/method.git",
        commit="xyz789",
        module_dir=Path(".repos/xyz789"),
        entrypoint=Path("run.py"),
        software_environment_id="python",
        resolved_environment=ResolvedEnvironment(
            backend_type=SoftwareBackendEnum.CONDA, reference="envs/python.yaml"
        ),
    )

    return [
        ResolvedNode(
            id="metrics-method1-default",
            stage_id="metrics",
            module_id="method1",
            param_id="default",
            module=module,
            outputs=["metrics/method1/default/scores.json"],
        ),
        ResolvedNode(
            id="metrics-method2-default",
            stage_id="metrics",
            module_id="method2",
            param_id="default",
            module=module,
            outputs=["metrics/method2/default/scores.json"],
        ),
        ResolvedNode(
            id="results-analyzer-default",
            stage_id="results",
            module_id="analyzer",
            param_id="default",
            module=module,
            outputs=["results/analyzer/default/data.csv"],
        ),
    ]


class TestResolveMetricCollectors:
    """Test the main resolve_metric_collectors function."""

    def test_resolve_single_collector(
        self, simple_collector, sample_resolved_nodes, mock_benchmark, mock_resolver
    ):
        """Test resolving a single collector without parameters."""
        collectors = [simple_collector]

        result = resolve_metric_collectors(
            metric_collectors=collectors,
            resolved_nodes=sample_resolved_nodes,
            benchmark=mock_benchmark,
            resolver=mock_resolver,
            quiet=True,
        )

        assert len(result) == 1
        node = result[0]

        # Check node properties
        assert node.id == "collector_plotting_default"
        assert node.stage_id == "_collector_plotting"
        assert node.module_id == "plotting"
        assert node.param_id == "default"

        # Check inputs gathered from stage
        assert len(node.inputs) == 2  # Two nodes in metrics stage
        assert "input_0" in node.inputs
        assert "input_1" in node.inputs

        # Check output path resolution
        assert len(node.outputs) == 1
        assert node.outputs[0] == "plotting/report.html"

    def test_resolve_collector_with_parameters(
        self,
        collector_with_params,
        sample_resolved_nodes,
        mock_benchmark,
        mock_resolver,
    ):
        """Test resolving a collector with parameter expansion."""
        collectors = [collector_with_params]

        result = resolve_metric_collectors(
            metric_collectors=collectors,
            resolved_nodes=sample_resolved_nodes,
            benchmark=mock_benchmark,
            resolver=mock_resolver,
            quiet=True,
        )

        # Should create 2 nodes (one per parameter value)
        assert len(result) == 2

        # Check parameter IDs are different
        param_ids = {node.param_id for node in result}
        assert "default" not in param_ids
        assert len(param_ids) == 2

        # Check output paths include param_id
        for node in result:
            assert node.param_id in node.outputs[0]

    def test_empty_collectors_list(
        self, sample_resolved_nodes, mock_benchmark, mock_resolver
    ):
        """Test handling empty collectors list."""
        result = resolve_metric_collectors(
            metric_collectors=[],
            resolved_nodes=sample_resolved_nodes,
            benchmark=mock_benchmark,
            resolver=mock_resolver,
            quiet=True,
        )

        assert result == []

    def test_collector_with_no_matching_stage(
        self, simple_collector, sample_resolved_nodes, mock_resolver
    ):
        """Test collector referencing non-existent stage."""
        # Create benchmark that returns None for stage lookup
        benchmark = Mock()
        benchmark.get_stage_by_output = Mock(return_value=None)
        benchmark.get_name = Mock(return_value="Test")
        benchmark.get_version = Mock(return_value="1.0")
        benchmark.get_author = Mock(return_value="Author")

        result = resolve_metric_collectors(
            metric_collectors=[simple_collector],
            resolved_nodes=sample_resolved_nodes,
            benchmark=benchmark,
            resolver=mock_resolver,
            quiet=True,
        )

        # Should still create node but with no inputs
        assert len(result) == 1
        assert len(result[0].inputs) == 0


class TestGatherCollectorInputs:
    """Test the _gather_collector_inputs function."""

    def test_gather_from_single_stage(
        self, simple_collector, sample_resolved_nodes, mock_benchmark
    ):
        """Test gathering inputs from a single stage."""
        inputs = _gather_collector_inputs(
            collector=simple_collector,
            resolved_nodes=sample_resolved_nodes,
            benchmark=mock_benchmark,
        )

        # Should gather from 2 nodes in metrics stage
        assert len(inputs) == 2
        assert "metrics/method1/default/scores.json" in inputs
        assert "metrics/method2/default/scores.json" in inputs

    def test_gather_from_multiple_stages(self, mock_benchmark, sample_resolved_nodes):
        """Test gathering inputs from multiple stages."""
        collector = MetricCollector(
            id="multi",
            name="Multi-stage collector",
            repository=Repository(url="http://test.git", commit="abc"),
            software_environment="python",
            inputs=["metrics.scores", "results.data"],
            outputs=[IOFile(id="out", path="{name}/out.html")],
        )

        inputs = _gather_collector_inputs(
            collector=collector,
            resolved_nodes=sample_resolved_nodes,
            benchmark=mock_benchmark,
        )

        # Should gather from both stages
        assert len(inputs) == 3
        assert "metrics/method1/default/scores.json" in inputs
        assert "metrics/method2/default/scores.json" in inputs
        assert "results/analyzer/default/data.csv" in inputs

    def test_gather_from_nonexistent_stage(
        self, simple_collector, sample_resolved_nodes
    ):
        """Test gathering when stage doesn't exist."""
        benchmark = Mock()
        benchmark.get_stage_by_output = Mock(return_value=None)

        inputs = _gather_collector_inputs(
            collector=simple_collector,
            resolved_nodes=sample_resolved_nodes,
            benchmark=benchmark,
        )

        assert inputs == []


class TestResolveCollectorOutputs:
    """Test the _resolve_collector_outputs function."""

    def test_simple_path_resolution(self):
        """Test basic path resolution with {name} substitution."""
        collector = MetricCollector(
            id="test_collector",
            name="Test",
            repository=Repository(url="http://test.git", commit="abc"),
            software_environment="python",
            inputs=["data"],
            outputs=[IOFile(id="report", path="{name}/report.html")],
        )

        outputs = _resolve_collector_outputs(collector, "default")

        assert len(outputs) == 1
        assert outputs[0] == "test_collector/report.html"

    def test_strip_input_placeholder(self):
        """Test stripping {input} placeholder."""
        collector = MetricCollector(
            id="plotter",
            name="Plotter",
            repository=Repository(url="http://test.git", commit="abc"),
            software_environment="python",
            inputs=["data"],
            outputs=[IOFile(id="plot", path="{input}/{name}/plot.png")],
        )

        outputs = _resolve_collector_outputs(collector, "default")

        assert outputs[0] == "plotter/plot.png"

    def test_path_with_parameters(self):
        """Test path resolution with parameter hash."""
        collector = MetricCollector(
            id="analyzer",
            name="Analyzer",
            repository=Repository(url="http://test.git", commit="abc"),
            software_environment="python",
            inputs=["data"],
            outputs=[IOFile(id="result", path="{name}/result.json")],
        )

        outputs = _resolve_collector_outputs(collector, "abc123")

        # Should insert param_id before filename
        assert outputs[0] == "analyzer/abc123/result.json"

    def test_clean_double_slashes(self):
        """Test cleaning up double slashes from path."""
        collector = MetricCollector(
            id="test",
            name="Test",
            repository=Repository(url="http://test.git", commit="abc"),
            software_environment="python",
            inputs=["data"],
            outputs=[IOFile(id="out", path="{input}//{name}//output.txt")],
        )

        outputs = _resolve_collector_outputs(collector, "default")

        # Should clean up double slashes
        assert "//" not in outputs[0]
        assert outputs[0] == "test/output.txt"


class TestExpandCollectorParameters:
    """Test the _expand_collector_parameters function."""

    def test_no_parameters(self, simple_collector):
        """Test expansion with no parameters."""
        result = _expand_collector_parameters(simple_collector)

        assert len(result) == 1
        assert result[0] is None

    def test_single_parameter(self):
        """Test expansion with single parameter."""
        collector = MetricCollector(
            id="test",
            name="Test",
            repository=Repository(url="http://test.git", commit="abc"),
            software_environment="python",
            inputs=["data"],
            outputs=[IOFile(id="out", path="out.txt")],
            parameters=[Parameter(id="alpha", values=[0.1, 0.5, 0.9])],
        )

        result = _expand_collector_parameters(collector)

        # Should expand to 3 parameter combinations
        assert len(result) == 3
        assert all(isinstance(p, Params) for p in result)

    def test_multiple_parameters(self):
        """Test expansion with multiple parameters (cartesian product)."""
        collector = MetricCollector(
            id="test",
            name="Test",
            repository=Repository(url="http://test.git", commit="abc"),
            software_environment="python",
            inputs=["data"],
            outputs=[IOFile(id="out", path="out.txt")],
            parameters=[
                Parameter(id="alpha", values=[0.1, 0.5]),
                Parameter(id="beta", values=[1, 2]),
            ],
        )

        result = _expand_collector_parameters(collector)

        # Should expand to 2 * 2 = 4 combinations
        assert len(result) == 4
        assert all(isinstance(p, Params) for p in result)
