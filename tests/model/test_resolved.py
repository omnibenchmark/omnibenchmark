"""
Tests for resolved entity classes (ResolvedNode, ResolvedModule, etc.)

These tests verify that the resolved entities correctly represent fully
resolved benchmark nodes with concrete paths and resolved parameters.
"""

import pytest
from pathlib import Path

from omnibenchmark.model import (
    ResolvedModule,
    ResolvedNode,
    ResolvedMetricCollector,
)
from omnibenchmark.benchmark.params import Params


class TestResolvedModule:
    """Tests for ResolvedModule class."""

    def test_create_resolved_module_with_absolute_paths(self):
        """Test creating a resolved module with absolute paths."""
        module = ResolvedModule(
            repository_url="https://github.com/test/module.git",
            commit="abc123",
            module_dir=Path("/home/user/.omnibenchmark/git/abc123"),
            entrypoint=Path("/home/user/.omnibenchmark/git/abc123/run.py"),
            software_environment_id="conda_env",
        )

        assert module.repository_url == "https://github.com/test/module.git"
        assert module.commit == "abc123"
        assert module.module_dir == Path("/home/user/.omnibenchmark/git/abc123")
        assert module.entrypoint == Path("/home/user/.omnibenchmark/git/abc123/run.py")
        assert module.software_environment_id == "conda_env"
        assert module.dirty is False

    def test_create_resolved_module_with_relative_paths(self):
        """Test that ResolvedModule allows relative paths."""
        # Relative paths are allowed (e.g., for cloning to cache)
        module = ResolvedModule(
            repository_url="https://github.com/test/module.git",
            commit="abc123",
            module_dir=Path(".omnibenchmark/git/abc123"),
            entrypoint=Path("run.py"),
            software_environment_id="conda_env",
        )

        assert module.module_dir == Path(".omnibenchmark/git/abc123")
        assert module.entrypoint == Path("run.py")

    def test_resolved_module_production_requires_commit(self):
        """Test that production mode (dirty=False) requires a commit."""
        with pytest.raises(
            ValueError, match="commit must be specified when dirty=False"
        ):
            ResolvedModule(
                repository_url="https://github.com/test/module.git",
                commit="",  # Empty commit
                module_dir=Path("/home/user/module"),
                entrypoint=Path("run.py"),
                software_environment_id="conda_env",
                dirty=False,  # Production mode
            )

    def test_resolved_module_dirty_allows_no_commit(self):
        """Test that dirty mode allows working copy without commit."""
        # In development, we can use working copies without specific commits
        module = ResolvedModule(
            repository_url="https://github.com/test/module.git",
            commit="",  # No commit (working copy)
            module_dir=Path("/home/user/dev/module"),
            entrypoint=Path("run.py"),
            software_environment_id="conda_env",
            dirty=True,  # Development mode
        )

        assert module.commit == ""
        assert module.dirty is True

    def test_resolved_module_is_immutable(self):
        """Test that ResolvedModule is immutable (frozen dataclass)."""
        module = ResolvedModule(
            repository_url="https://github.com/test/module.git",
            commit="abc123",
            module_dir=Path("/home/user/.omnibenchmark/git/abc123"),
            entrypoint=Path("/home/user/.omnibenchmark/git/abc123/run.py"),
            software_environment_id="conda_env",
        )

        with pytest.raises(Exception):  # FrozenInstanceError
            module.commit = "def456"


class TestResolvedNode:
    """Tests for ResolvedNode class."""

    def test_create_simple_resolved_node(self):
        """Test creating a simple resolved node without parameters."""
        module = ResolvedModule(
            repository_url="https://github.com/test/module.git",
            commit="abc123",
            module_dir=Path("/home/user/.omnibenchmark/git/abc123"),
            entrypoint=Path("/home/user/.omnibenchmark/git/abc123/run.py"),
            software_environment_id="conda_env",
        )

        node = ResolvedNode(
            id="data-D1-default",
            stage_id="data",
            module_id="D1",
            param_id="default",
            module=module,
            outputs=["{input}/data/D1/default/{dataset}.json"],
        )

        assert node.id == "data-D1-default"
        assert node.stage_id == "data"
        assert node.module_id == "D1"
        assert node.param_id == "default"
        assert node.module == module
        assert node.outputs == ["{input}/data/D1/default/{dataset}.json"]
        assert node.is_entrypoint()  # No inputs

    def test_resolved_node_with_parameters(self):
        """Test resolved node with parameters."""
        module = ResolvedModule(
            repository_url="https://github.com/test/module.git",
            commit="abc123",
            module_dir=Path("/home/user/.omnibenchmark/git/abc123"),
            entrypoint=Path("/home/user/.omnibenchmark/git/abc123/run.py"),
            software_environment_id="conda_env",
        )

        params = Params({"method": "cosine", "threshold": 0.1})

        node = ResolvedNode(
            id="methods-M1-.12345678",
            stage_id="methods",
            module_id="M1",
            param_id=".12345678",
            module=module,
            parameters=params,
            param_dir_template="{input}/methods/M1/.12345678",
            param_symlink_template="{input}/methods/M1/method-cosine_threshold-0.1",
            outputs=["{input}/methods/M1/.12345678/{dataset}.result.json"],
        )

        assert node.parameters == params
        assert node.param_dir_template == "{input}/methods/M1/.12345678"
        assert (
            node.param_symlink_template
            == "{input}/methods/M1/method-cosine_threshold-0.1"
        )

        # Test parameter methods
        assert node.get_parameter_hash() == params.hash_short()
        assert node.get_parameter_json() == params.serialize()
        assert "--method" in node.get_parameter_cli_args()

    def test_resolved_node_with_inputs_and_parent(self):
        """Test resolved node with inputs and parent reference."""
        module = ResolvedModule(
            repository_url="https://github.com/test/module.git",
            commit="abc123",
            module_dir=Path("/home/user/.omnibenchmark/git/abc123"),
            entrypoint=Path("/home/user/.omnibenchmark/git/abc123/run.py"),
            software_environment_id="conda_env",
        )

        node = ResolvedNode(
            id="methods-M1-default-after_data",
            stage_id="methods",
            module_id="M1",
            param_id="default",
            module=module,
            parent_id="data-D1-default",
            inputs={
                "data.raw": "{input}/data/D1/default/{dataset}.json",
            },
            outputs=["{input}/methods/M1/default/{dataset}.result.json"],
        )

        assert node.parent_id == "data-D1-default"
        assert node.inputs == {"data.raw": "{input}/data/D1/default/{dataset}.json"}
        assert not node.is_entrypoint()  # Has inputs

        # Test input accessors
        assert node.get_input_dict() == {
            "data.raw": "{input}/data/D1/default/{dataset}.json"
        }
        assert node.get_input_list() == ["{input}/data/D1/default/{dataset}.json"]

    def test_resolved_node_id_creation(self):
        """Test node ID creation utility."""
        assert ResolvedNode.create_id("data", "D1", "default") == "data-D1-default"
        assert (
            ResolvedNode.create_id("methods", "M1", ".abc123", "data")
            == "methods-M1-.abc123-after_data"
        )

    def test_resolved_node_to_dict(self):
        """Test conversion to dictionary for serialization."""
        module = ResolvedModule(
            repository_url="https://github.com/test/module.git",
            commit="abc123",
            module_dir=Path("/home/user/.omnibenchmark/git/abc123"),
            entrypoint=Path("/home/user/.omnibenchmark/git/abc123/run.py"),
            software_environment_id="conda_env",
        )

        params = Params({"method": "test"})

        node = ResolvedNode(
            id="stage-module-.12345678",
            stage_id="stage",
            module_id="module",
            param_id=".12345678",
            module=module,
            parameters=params,
            benchmark_name="TestBench",
            benchmark_version="1.0",
            benchmark_author="Tester",
        )

        d = node.to_dict()

        assert d["id"] == "stage-module-.12345678"
        assert d["stage_id"] == "stage"
        assert d["module_id"] == "module"
        assert d["param_id"] == ".12345678"
        assert d["module"]["repository_url"] == "https://github.com/test/module.git"
        assert d["module"]["commit"] == "abc123"
        assert d["benchmark_name"] == "TestBench"
        assert d["benchmark_version"] == "1.0"
        assert d["benchmark_author"] == "Tester"

    def test_resolved_node_is_immutable(self):
        """Test that ResolvedNode is immutable."""
        module = ResolvedModule(
            repository_url="https://github.com/test/module.git",
            commit="abc123",
            module_dir=Path("/home/user/.omnibenchmark/git/abc123"),
            entrypoint=Path("/home/user/.omnibenchmark/git/abc123/run.py"),
            software_environment_id="conda_env",
        )

        node = ResolvedNode(
            id="data-D1-default",
            stage_id="data",
            module_id="D1",
            param_id="default",
            module=module,
        )

        with pytest.raises(Exception):  # FrozenInstanceError
            node.stage_id = "other_stage"

    def test_resolved_node_without_parameters(self):
        """Test node without parameters returns sensible defaults."""
        module = ResolvedModule(
            repository_url="https://github.com/test/module.git",
            commit="abc123",
            module_dir=Path("/home/user/.omnibenchmark/git/abc123"),
            entrypoint=Path("/home/user/.omnibenchmark/git/abc123/run.py"),
            software_environment_id="conda_env",
        )

        node = ResolvedNode(
            id="data-D1-default",
            stage_id="data",
            module_id="D1",
            param_id="default",
            module=module,
        )

        assert node.get_parameter_cli_args() == []
        assert node.get_parameter_json() == "{}"
        assert node.get_parameter_hash() == "default"

    def test_resolved_node_string_representation(self):
        """Test string representations."""
        module = ResolvedModule(
            repository_url="https://github.com/test/module.git",
            commit="abc123",
            module_dir=Path("/home/user/.omnibenchmark/git/abc123"),
            entrypoint=Path("/home/user/.omnibenchmark/git/abc123/run.py"),
            software_environment_id="conda_env",
        )

        node = ResolvedNode(
            id="data-D1-default",
            stage_id="data",
            module_id="D1",
            param_id="default",
            module=module,
        )

        assert str(node) == "data-D1-default"
        assert repr(node) == "ResolvedNode(data-D1-default)"


class TestResolvedMetricCollector:
    """Tests for ResolvedMetricCollector class."""

    def test_create_resolved_metric_collector(self):
        """Test creating a resolved metric collector."""
        module = ResolvedModule(
            repository_url="https://github.com/test/collector.git",
            commit="xyz789",
            module_dir=Path("/home/user/.omnibenchmark/git/xyz789"),
            entrypoint=Path("/home/user/.omnibenchmark/git/xyz789/collect.py"),
            software_environment_id="conda_env",
        )

        collector = ResolvedMetricCollector(
            id="MC1",
            name="Test Metric Collector",
            module=module,
            input_patterns=["{input}/methods/*/default/*.result.json"],
            outputs=["{input}/metrics/metrics.json"],
        )

        assert collector.id == "MC1"
        assert collector.name == "Test Metric Collector"
        assert collector.module == module
        assert collector.input_patterns == ["{input}/methods/*/default/*.result.json"]
        assert collector.outputs == ["{input}/metrics/metrics.json"]

    def test_metric_collector_accessors(self):
        """Test accessor methods for metric collector."""
        module = ResolvedModule(
            repository_url="https://github.com/test/collector.git",
            commit="xyz789",
            module_dir=Path("/home/user/.omnibenchmark/git/xyz789"),
            entrypoint=Path("/home/user/.omnibenchmark/git/xyz789/collect.py"),
            software_environment_id="conda_env",
        )

        collector = ResolvedMetricCollector(
            id="MC1",
            module=module,
            input_patterns=["pattern1", "pattern2"],
            outputs=["output1", "output2"],
        )

        assert collector.get_input_list() == ["pattern1", "pattern2"]
        assert collector.get_output_list() == ["output1", "output2"]

    def test_metric_collector_to_dict(self):
        """Test conversion to dictionary."""
        module = ResolvedModule(
            repository_url="https://github.com/test/collector.git",
            commit="xyz789",
            module_dir=Path("/home/user/.omnibenchmark/git/xyz789"),
            entrypoint=Path("/home/user/.omnibenchmark/git/xyz789/collect.py"),
            software_environment_id="conda_env",
        )

        collector = ResolvedMetricCollector(
            id="MC1",
            name="Test Collector",
            module=module,
            benchmark_name="TestBench",
        )

        d = collector.to_dict()

        assert d["id"] == "MC1"
        assert d["name"] == "Test Collector"
        assert d["module"]["repository_url"] == "https://github.com/test/collector.git"
        assert d["benchmark_name"] == "TestBench"

    def test_metric_collector_string_representation(self):
        """Test string representations."""
        module = ResolvedModule(
            repository_url="https://github.com/test/collector.git",
            commit="xyz789",
            module_dir=Path("/home/user/.omnibenchmark/git/xyz789"),
            entrypoint=Path("/home/user/.omnibenchmark/git/xyz789/collect.py"),
            software_environment_id="conda_env",
        )

        collector = ResolvedMetricCollector(
            id="MC1",
            module=module,
        )

        assert str(collector) == "MC1"
        assert repr(collector) == "ResolvedMetricCollector(MC1)"


class TestResolvedNodeParameterHandling:
    """Tests for parameter handling in ResolvedNode."""

    def test_parameter_cli_args_gnu_style(self):
        """Test generating GNU-style CLI arguments."""
        module = ResolvedModule(
            repository_url="https://github.com/test/module.git",
            commit="abc123",
            module_dir=Path("/home/user/.omnibenchmark/git/abc123"),
            entrypoint=Path("/home/user/.omnibenchmark/git/abc123/run.py"),
            software_environment_id="conda_env",
        )

        params = Params({"method": "cosine", "k": 10})

        node = ResolvedNode(
            id="test-node-default",
            stage_id="test",
            module_id="node",
            param_id="default",
            module=module,
            parameters=params,
        )

        args = node.get_parameter_cli_args(style="gnu")
        assert "--k" in args
        assert "10" in args
        assert "--method" in args
        assert "cosine" in args

    def test_parameter_cli_args_equals_style(self):
        """Test generating equals-style CLI arguments."""
        module = ResolvedModule(
            repository_url="https://github.com/test/module.git",
            commit="abc123",
            module_dir=Path("/home/user/.omnibenchmark/git/abc123"),
            entrypoint=Path("/home/user/.omnibenchmark/git/abc123/run.py"),
            software_environment_id="conda_env",
        )

        params = Params({"method": "cosine"})

        node = ResolvedNode(
            id="test-node-default",
            stage_id="test",
            module_id="node",
            param_id="default",
            module=module,
            parameters=params,
        )

        args = node.get_parameter_cli_args(style="equals")
        assert "--method=cosine" in args
