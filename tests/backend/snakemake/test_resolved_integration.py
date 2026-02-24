"""
Integration tests demonstrating ResolvedNode usage with e2e YAML examples.

These tests show how to build ResolvedNode entities from parsed benchmark
YAML files, simulating the resolution phase of the execution pipeline.
"""

import pytest
from pathlib import Path

from omnibenchmark.model import (
    Benchmark,
    ResolvedModule,
    ResolvedNode,
    ResolvedMetricCollector,
)
from omnibenchmark.benchmark.params import Params


@pytest.fixture
def simple_benchmark_yaml():
    """Path to simple e2e benchmark YAML."""
    return (
        Path(__file__).parent.parent.parent
        / "e2e"
        / "configs"
        / "01_data_and_methods.yaml"
    )


@pytest.fixture
def benchmark_with_metrics_yaml():
    """Path to benchmark with metrics YAML."""
    return (
        Path(__file__).parent.parent.parent
        / "e2e"
        / "configs"
        / "02_data_methods_metrics.yaml"
    )


class TestResolvedNodeFromBenchmark:
    """Tests showing how to create ResolvedNode from parsed Benchmark."""

    def test_load_simple_benchmark(self, simple_benchmark_yaml):
        """Test loading a simple benchmark YAML."""
        benchmark = Benchmark.from_yaml(simple_benchmark_yaml)

        assert benchmark.id == "LinearArithmeticDataset"
        assert len(benchmark.stages) == 2
        assert benchmark.stages[0].id == "data"
        assert benchmark.stages[1].id == "methods"

    def test_create_resolved_node_for_data_stage(self, simple_benchmark_yaml):
        """Test creating a ResolvedNode for a data stage module."""
        benchmark = Benchmark.from_yaml(simple_benchmark_yaml)

        # Get first stage (data)
        data_stage = benchmark.stages[0]
        assert data_stage.id == "data"

        # Get first module (D1)
        module = data_stage.modules[0]
        assert module.id == "D1"

        # Simulate module resolution (in real code this would clone the repo)
        resolved_module = ResolvedModule(
            repository_url=module.repository.url,
            commit=module.repository.commit,
            module_dir=Path("/home/user/.omnibenchmark/git/4ff8427"),
            entrypoint=Path("/home/user/.omnibenchmark/git/4ff8427/run.py"),
            software_environment_id=module.software_environment,
        )

        # Create ResolvedNode
        node = ResolvedNode(
            id="data-D1-default",
            stage_id=data_stage.id,
            module_id=module.id,
            param_id="default",
            module=resolved_module,
            outputs=["{input}/data/D1/default/{dataset}_data.json"],
            benchmark_name=benchmark.get_name(),
            benchmark_version=benchmark.get_version(),
            benchmark_author=benchmark.get_author(),
        )

        assert node.id == "data-D1-default"
        assert node.stage_id == "data"
        assert node.module_id == "D1"
        assert node.is_entrypoint()  # Data stage has no inputs
        assert len(node.outputs) == 1

    def test_create_resolved_node_with_parameters(self, simple_benchmark_yaml):
        """Test creating a ResolvedNode with parameters."""
        benchmark = Benchmark.from_yaml(simple_benchmark_yaml)

        # Get data stage module D1 which has parameters
        data_stage = benchmark.stages[0]
        module = data_stage.modules[0]

        # Get parameters from module
        assert module.parameters is not None
        assert len(module.parameters) > 0

        param = module.parameters[0]
        # This is the old format with values: ["--evaluate", "1+1"]
        params_list = Params.expand_from_parameter(param)
        assert len(params_list) > 0

        params = params_list[0]

        # Simulate module resolution
        resolved_module = ResolvedModule(
            repository_url=module.repository.url,
            commit=module.repository.commit,
            module_dir=Path("/home/user/.omnibenchmark/git/4ff8427"),
            entrypoint=Path("/home/user/.omnibenchmark/git/4ff8427/run.py"),
            software_environment_id=module.software_environment,
        )

        # Create ResolvedNode with parameters
        param_hash = params.hash_short()
        node = ResolvedNode(
            id=f"data-D1-.{param_hash}",
            stage_id=data_stage.id,
            module_id=module.id,
            param_id=f".{param_hash}",
            module=resolved_module,
            parameters=params,
            param_dir_template=f"{{input}}/data/D1/.{param_hash}",
            outputs=[f"{{input}}/data/D1/.{param_hash}/{{dataset}}_data.json"],
            benchmark_name=benchmark.get_name(),
            benchmark_version=benchmark.get_version(),
            benchmark_author=benchmark.get_author(),
        )

        assert node.parameters is not None
        assert node.get_parameter_hash() == param_hash
        assert node.param_dir_template.endswith(param_hash)

    def test_create_resolved_node_for_methods_stage(self, simple_benchmark_yaml):
        """Test creating a ResolvedNode for a methods stage with inputs."""
        benchmark = Benchmark.from_yaml(simple_benchmark_yaml)

        # Get methods stage
        methods_stage = benchmark.stages[1]
        assert methods_stage.id == "methods"

        # Get first module (M1)
        module = methods_stage.modules[0]
        assert module.id == "M1"

        # Methods stage has inputs from data stage
        assert methods_stage.inputs is not None
        assert len(methods_stage.inputs) > 0

        # Simulate module resolution
        resolved_module = ResolvedModule(
            repository_url=module.repository.url,
            commit=module.repository.commit,
            module_dir=Path("/home/user/.omnibenchmark/git/4ff8427"),
            entrypoint=Path("/home/user/.omnibenchmark/git/4ff8427/run.py"),
            software_environment_id=module.software_environment,
        )

        # Create ResolvedNode with inputs
        node = ResolvedNode(
            id="methods-M1-default-after_data",
            stage_id=methods_stage.id,
            module_id=module.id,
            param_id="default",
            module=resolved_module,
            parent_id="data-D1-default",
            inputs={
                "data.raw": "{input}/data/D1/default/{dataset}_data.json",
            },
            outputs=["{input}/methods/M1/default/{dataset}_data.json"],
            benchmark_name=benchmark.get_name(),
            benchmark_version=benchmark.get_version(),
            benchmark_author=benchmark.get_author(),
        )

        assert node.id == "methods-M1-default-after_data"
        assert not node.is_entrypoint()  # Has inputs
        assert node.parent_id == "data-D1-default"
        assert "data.raw" in node.inputs
        assert node.get_input_list() == ["{input}/data/D1/default/{dataset}_data.json"]


class TestResolvedMetricCollectorFromBenchmark:
    """Tests showing how to create ResolvedMetricCollector from parsed Benchmark."""

    def test_create_resolved_metric_collector(self, benchmark_with_metrics_yaml):
        """Test creating a ResolvedMetricCollector from benchmark with metrics."""
        benchmark = Benchmark.from_yaml(benchmark_with_metrics_yaml)

        # Check that benchmark has metric collectors
        assert benchmark.metric_collectors is not None
        assert len(benchmark.metric_collectors) > 0

        # Get first metric collector
        mc = benchmark.metric_collectors[0]
        assert mc.id == "MC1"

        # Simulate module resolution for metric collector
        resolved_module = ResolvedModule(
            repository_url=mc.repository.url,
            commit=mc.repository.commit,
            module_dir=Path("/home/user/.omnibenchmark/git/9e3a99a"),
            entrypoint=Path("/home/user/.omnibenchmark/git/9e3a99a/collect.py"),
            software_environment_id=mc.software_environment,
        )

        # Create ResolvedMetricCollector
        collector = ResolvedMetricCollector(
            id=mc.id,
            name=mc.name,
            module=resolved_module,
            input_patterns=["{input}/methods/*/default/*.json"],
            outputs=["{input}/metrics/metrics.json"],
            benchmark_name=benchmark.get_name(),
            benchmark_version=benchmark.get_version(),
            benchmark_author=benchmark.get_author(),
        )

        assert collector.id == "MC1"
        assert collector.name == "Gather and calculate metrics from method results"
        assert len(collector.input_patterns) > 0
        assert len(collector.outputs) > 0


class TestResolvedNodeSerialization:
    """Tests for serialization and conversion of resolved nodes."""

    def test_resolved_node_to_dict_roundtrip(self, simple_benchmark_yaml):
        """Test converting ResolvedNode to dict and back."""
        benchmark = Benchmark.from_yaml(simple_benchmark_yaml)
        data_stage = benchmark.stages[0]
        module = data_stage.modules[0]

        # Create a resolved module
        resolved_module = ResolvedModule(
            repository_url=module.repository.url,
            commit=module.repository.commit,
            module_dir=Path("/home/user/.omnibenchmark/git/4ff8427"),
            entrypoint=Path("/home/user/.omnibenchmark/git/4ff8427/run.py"),
            software_environment_id=module.software_environment,
        )

        # Create a resolved node
        node = ResolvedNode(
            id="data-D1-default",
            stage_id=data_stage.id,
            module_id=module.id,
            param_id="default",
            module=resolved_module,
            outputs=["{input}/data/D1/default/{dataset}_data.json"],
            benchmark_name=benchmark.get_name(),
            benchmark_version=benchmark.get_version(),
            benchmark_author=benchmark.get_author(),
        )

        # Convert to dict
        node_dict = node.to_dict()

        # Verify dict structure
        assert node_dict["id"] == "data-D1-default"
        assert node_dict["stage_id"] == "data"
        assert node_dict["module_id"] == "D1"
        assert "module" in node_dict
        assert node_dict["module"]["repository_url"] == module.repository.url
        assert node_dict["module"]["commit"] == module.repository.commit
        assert node_dict["benchmark_name"] == benchmark.get_name()

        # Verify it's JSON-serializable (all values are primitive types)
        import json

        json_str = json.dumps(node_dict)
        assert len(json_str) > 0


class TestResolvedNodeDesignPrinciples:
    """Tests verifying the design principles of resolved entities."""

    def test_resolved_entities_are_immutable(self):
        """Verify that resolved entities are immutable (frozen dataclasses)."""
        module = ResolvedModule(
            repository_url="https://github.com/test/module.git",
            commit="abc123",
            module_dir=Path("/home/user/.omnibenchmark/git/abc123"),
            entrypoint=Path("/home/user/.omnibenchmark/git/abc123/run.py"),
            software_environment_id="host",
        )

        node = ResolvedNode(
            id="test-node-default",
            stage_id="test",
            module_id="node",
            param_id="default",
            module=module,
        )

        # Verify immutability
        with pytest.raises(Exception):
            node.stage_id = "different"

        with pytest.raises(Exception):
            module.commit = "different"

    def test_resolved_node_no_backend_specific_methods(self):
        """Verify ResolvedNode has no backend-specific methods."""
        # Get all public methods
        methods = [m for m in dir(ResolvedNode) if not m.startswith("_")]

        # These are all generic data access methods, not backend-specific
        backend_keywords = ["snakemake", "snake", "workflow", "rule"]
        for method in methods:
            for keyword in backend_keywords:
                assert (
                    keyword not in method.lower()
                ), f"Found backend-specific method: {method}"

    def test_resolved_module_production_requires_commit(self):
        """Verify ResolvedModule enforces commit in production mode."""
        # Production mode (dirty=False) requires a commit
        with pytest.raises(
            ValueError, match="commit must be specified when dirty=False"
        ):
            ResolvedModule(
                repository_url="https://github.com/test/module.git",
                commit="",  # Empty commit
                module_dir=Path("relative/path"),
                entrypoint=Path("run.py"),
                software_environment_id="host",
                dirty=False,  # Production mode
            )

    def test_resolved_module_allows_relative_paths(self):
        """Verify ResolvedModule allows relative paths for cache."""
        # Relative paths are allowed (used for git cache)
        module = ResolvedModule(
            repository_url="https://github.com/test/module.git",
            commit="abc123",
            module_dir=Path(".omnibenchmark/git/abc123"),
            entrypoint=Path("run.py"),
            software_environment_id="host",
        )
        assert module.module_dir == Path(".omnibenchmark/git/abc123")
        assert module.entrypoint == Path("run.py")
