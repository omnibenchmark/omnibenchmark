"""
Unit tests for DebugSnakemakeGenerator (backend/_snakemake_debug.py).

Tests the dry-run Snakefile generator that emits echo/touch commands instead
of real execution logic for testing DAG structure and path resolution.

All tests marked as 'short' since they test pure logic without external dependencies.
"""

import pytest
from io import StringIO
from pathlib import Path

from omnibenchmark.backend._snakemake_debug import DebugSnakemakeGenerator
from omnibenchmark.model.resolved import (
    ResolvedNode,
    ResolvedModule,
    ResolvedEnvironment,
)
from omnibenchmark.model import SoftwareBackendEnum
from omnibenchmark.benchmark.params import Params


@pytest.fixture
def basic_resolved_module():
    """Create a basic ResolvedModule for testing."""
    return ResolvedModule(
        repository_url="https://github.com/test/module.git",
        commit="abc123",
        module_dir=Path(".repos/abc123"),
        entrypoint=Path("run.py"),
        software_environment_id="python",
        resolved_environment=ResolvedEnvironment(
            backend_type=SoftwareBackendEnum.conda,
            reference="envs/python.yaml",
        ),
        has_shebang=False,
        interpreter="python3",
    )


@pytest.fixture
def basic_node(basic_resolved_module):
    """Create a basic ResolvedNode without inputs or parameters."""
    return ResolvedNode(
        id="test_node_default",
        stage_id="test_stage",
        module_id="test_module",
        param_id="default",
        module=basic_resolved_module,
        outputs=["output/result.txt"],
    )


@pytest.fixture
def node_with_inputs(basic_resolved_module):
    """Create a ResolvedNode with inputs."""
    return ResolvedNode(
        id="processor_default",
        stage_id="processing",
        module_id="processor",
        param_id="default",
        module=basic_resolved_module,
        inputs={"data": "input/data.csv", "config": "input/config.json"},
        outputs=["output/processed.csv"],
        input_name_mapping={"data": "data", "config": "config"},
    )


@pytest.fixture
def node_with_parameters(basic_resolved_module):
    """Create a ResolvedNode with parameters."""
    params = Params({"alpha": 0.05, "beta": 1.5})
    return ResolvedNode(
        id="analyzer_abc123",
        stage_id="analysis",
        module_id="analyzer",
        param_id="abc123",
        module=basic_resolved_module,
        outputs=["output/analysis.txt"],
        parameters=params,
    )


@pytest.fixture
def gather_node(basic_resolved_module):
    """Create a gather/collector node."""
    return ResolvedNode(
        id="collector_plotting_default",
        stage_id="_collector_plotting",
        module_id="plotting",
        param_id="default",
        module=basic_resolved_module,
        inputs={
            "input_0": "metrics/method1/scores.json",
            "input_1": "metrics/method2/scores.json",
            "input_2": "metrics/method3/scores.json",
        },
        outputs=["plotting/report.html"],
        input_name_mapping={
            "input_0": "scores",
            "input_1": "scores",
            "input_2": "scores",
        },
        is_gather=True,
    )


@pytest.fixture
def collector_node(basic_resolved_module):
    """Create a metric collector node (is_collector=True)."""
    return ResolvedNode(
        id="collector_metrics_default",
        stage_id="_collector_metrics",
        module_id="metrics",
        param_id="default",
        module=basic_resolved_module,
        inputs={
            "input_0": "results/a.csv",
            "input_1": "results/b.csv",
        },
        outputs=["metrics/summary.json"],
        input_name_mapping={
            "input_0": "results",
            "input_1": "results",
        },
        is_collector=True,
    )


@pytest.mark.short
class TestDebugSnakemakeGenerator:
    """Test DebugSnakemakeGenerator instantiation and basic properties."""

    def test_instantiation(self):
        """Test that DebugSnakemakeGenerator can be instantiated."""
        generator = DebugSnakemakeGenerator(
            benchmark_name="test",
            benchmark_version="1.0",
            benchmark_author="author",
        )
        assert generator is not None
        assert isinstance(generator, DebugSnakemakeGenerator)

    def test_inherits_from_snakemake_generator(self):
        """Test that DebugSnakemakeGenerator inherits from SnakemakeGenerator."""
        from omnibenchmark.backend.snakemake import SnakemakeGenerator

        generator = DebugSnakemakeGenerator(
            benchmark_name="test",
            benchmark_version="1.0",
            benchmark_author="author",
        )
        assert isinstance(generator, SnakemakeGenerator)


@pytest.mark.short
class TestWriteShellBasicNode:
    """Test _write_shell method for basic nodes."""

    def test_basic_node_without_inputs_or_params(self, basic_node):
        """Test shell block for node with only outputs."""
        generator = DebugSnakemakeGenerator(
            benchmark_name="test",
            benchmark_version="1.0",
            benchmark_author="author",
        )

        output = StringIO()
        generator._write_shell(output, basic_node)
        content = output.getvalue()

        # Should contain rule identification
        assert f'echo "RULE: {basic_node.id}"' in content
        assert 'echo "Module: {params.module_dir}"' in content
        assert 'echo "Entrypoint: {params.entrypoint}"' in content

        # Should contain output section
        assert 'echo "Outputs:"' in content
        assert 'echo "  {output[0]}"' in content

        # Should contain execution preview
        assert 'echo "Would execute:"' in content
        assert f'echo "    --name {basic_node.module_id} \\"' in content

        # Should create output directory and touch files
        assert "mkdir -p $(dirname {output[0]})" in content
        assert "touch {output[0]}" in content

    def test_node_with_inputs(self, node_with_inputs):
        """Test shell block includes input echoing."""
        generator = DebugSnakemakeGenerator(
            benchmark_name="test",
            benchmark_version="1.0",
            benchmark_author="author",
        )

        output = StringIO()
        generator._write_shell(output, node_with_inputs)
        content = output.getvalue()

        # Should contain input section
        assert 'echo "Inputs:"' in content
        assert 'echo "  data: {input.data}"' in content
        assert 'echo "  config: {input.config}"' in content

        # Should map inputs back to original names
        assert 'echo "    --data {input.data} \\"' in content
        assert 'echo "    --config {input.config} \\"' in content

    def test_node_with_parameters(self, node_with_parameters):
        """Test shell block includes parameter handling."""
        generator = DebugSnakemakeGenerator(
            benchmark_name="test",
            benchmark_version="1.0",
            benchmark_author="author",
        )

        output = StringIO()
        generator._write_shell(output, node_with_parameters)
        content = output.getvalue()

        # Should contain parameters section
        assert 'echo "Parameters: {params.cli_args}"' in content

        # Should include parameters in execution preview
        assert 'echo "    {params.cli_args}"' in content

        # Should create parameter symlink
        assert "ln -sfn" in content
        assert node_with_parameters.parameters.hash_short() in content

    def test_node_with_multiple_outputs(self, basic_resolved_module):
        """Test shell block handles multiple outputs."""
        node = ResolvedNode(
            id="multi_output",
            stage_id="stage",
            module_id="module",
            param_id="default",
            module=basic_resolved_module,
            outputs=["output/file1.txt", "output/file2.csv", "output/file3.json"],
        )

        generator = DebugSnakemakeGenerator(
            benchmark_name="test",
            benchmark_version="1.0",
            benchmark_author="author",
        )

        output = StringIO()
        generator._write_shell(output, node)
        content = output.getvalue()

        # Should echo all outputs
        assert 'echo "  {output[0]}"' in content
        assert 'echo "  {output[1]}"' in content
        assert 'echo "  {output[2]}"' in content

        # Should touch all outputs
        assert "touch {output[0]}" in content
        assert "touch {output[1]}" in content
        assert "touch {output[2]}" in content

    def test_empty_inputs_dict(self, basic_resolved_module):
        """Test node with empty inputs dict doesn't crash."""
        node = ResolvedNode(
            id="no_inputs",
            stage_id="stage",
            module_id="module",
            param_id="default",
            module=basic_resolved_module,
            inputs={},
            outputs=["output/result.txt"],
        )

        generator = DebugSnakemakeGenerator(
            benchmark_name="test",
            benchmark_version="1.0",
            benchmark_author="author",
        )

        output = StringIO()
        generator._write_shell(output, node)
        content = output.getvalue()

        # Should not have inputs section
        assert 'echo "Inputs:"' not in content


@pytest.mark.short
class TestWriteGatherShell:
    """Test _write_gather_shell method for gather/collector nodes."""

    def test_gather_node_basic(self, gather_node):
        """Test gather shell block for basic gather node."""
        generator = DebugSnakemakeGenerator(
            benchmark_name="test",
            benchmark_version="1.0",
            benchmark_author="author",
        )

        output = StringIO()
        generator._write_gather_shell(output, gather_node)
        content = output.getvalue()

        # Should identify as gather stage
        assert 'echo "GATHER STAGE: plotting"' in content
        assert 'echo "Module: {params.module_dir}"' in content
        assert 'echo "Entrypoint: {params.entrypoint}"' in content

    def test_collector_node_identification(self, collector_node):
        """Test that collector nodes are identified correctly."""
        generator = DebugSnakemakeGenerator(
            benchmark_name="test",
            benchmark_version="1.0",
            benchmark_author="author",
        )

        output = StringIO()
        generator._write_gather_shell(output, collector_node)
        content = output.getvalue()

        # Should identify as metric collector
        assert 'echo "METRIC COLLECTOR: metrics"' in content

    def test_gather_inputs_grouped_by_name(self, gather_node):
        """Test that gather inputs are grouped by their original names."""
        generator = DebugSnakemakeGenerator(
            benchmark_name="test",
            benchmark_version="1.0",
            benchmark_author="author",
        )

        output = StringIO()
        generator._write_gather_shell(output, gather_node)
        content = output.getvalue()

        # Should have inputs by name section
        assert 'echo "Inputs (by name):"' in content
        assert 'echo "  --scores:"' in content

        # Should list all input keys under the name
        assert 'echo "    {input.input_0}"' in content
        assert 'echo "    {input.input_1}"' in content
        assert 'echo "    {input.input_2}"' in content

    def test_gather_with_multiple_input_names(self, basic_resolved_module):
        """Test gather node with multiple distinct input names."""
        node = ResolvedNode(
            id="multi_gather",
            stage_id="_collector_multi",
            module_id="multi",
            param_id="default",
            module=basic_resolved_module,
            inputs={
                "input_0": "data/a.csv",
                "input_1": "data/b.csv",
                "input_2": "metrics/x.json",
                "input_3": "metrics/y.json",
            },
            outputs=["summary.html"],
            input_name_mapping={
                "input_0": "data",
                "input_1": "data",
                "input_2": "metrics",
                "input_3": "metrics",
            },
            is_gather=True,
        )

        generator = DebugSnakemakeGenerator(
            benchmark_name="test",
            benchmark_version="1.0",
            benchmark_author="author",
        )

        output = StringIO()
        generator._write_gather_shell(output, node)
        content = output.getvalue()

        # Should have both input names
        assert 'echo "  --data:"' in content
        assert 'echo "  --metrics:"' in content

        # Should show grouped inputs
        assert 'echo "    {input.input_0}"' in content
        assert 'echo "    {input.input_1}"' in content
        assert 'echo "    {input.input_2}"' in content
        assert 'echo "    {input.input_3}"' in content

    def test_gather_outputs_handling(self, gather_node):
        """Test that gather outputs are properly echoed and touched."""
        generator = DebugSnakemakeGenerator(
            benchmark_name="test",
            benchmark_version="1.0",
            benchmark_author="author",
        )

        output = StringIO()
        generator._write_gather_shell(output, gather_node)
        content = output.getvalue()

        # Should echo outputs
        assert 'echo "Outputs:"' in content
        assert 'echo "  {output[0]}"' in content

        # Should touch outputs
        assert "touch {output[0]}" in content

    def test_gather_with_parameters(self, basic_resolved_module):
        """Test gather node with parameters."""
        params = Params({"threshold": 0.05})
        node = ResolvedNode(
            id="gather_with_params",
            stage_id="_collector_test",
            module_id="test",
            param_id="abc123",
            module=basic_resolved_module,
            inputs={"input_0": "data.csv"},
            outputs=["result.txt"],
            input_name_mapping={"input_0": "data"},
            parameters=params,
            is_gather=True,
        )

        generator = DebugSnakemakeGenerator(
            benchmark_name="test",
            benchmark_version="1.0",
            benchmark_author="author",
        )

        output = StringIO()
        generator._write_gather_shell(output, node)
        content = output.getvalue()

        # Should include parameters
        assert 'echo "Parameters: {params.cli_args}"' in content
        assert 'echo "    {params.cli_args}"' in content

    def test_gather_without_inputs(self, basic_resolved_module):
        """Test gather node without inputs (edge case)."""
        node = ResolvedNode(
            id="empty_gather",
            stage_id="_collector_empty",
            module_id="empty",
            param_id="default",
            module=basic_resolved_module,
            inputs={},
            outputs=["output.txt"],
            is_gather=True,
        )

        generator = DebugSnakemakeGenerator(
            benchmark_name="test",
            benchmark_version="1.0",
            benchmark_author="author",
        )

        output = StringIO()
        generator._write_gather_shell(output, node)
        content = output.getvalue()

        # Should not crash and should still create outputs
        assert "touch {output[0]}" in content

    def test_gather_execution_preview(self, gather_node):
        """Test that gather nodes show proper execution preview."""
        generator = DebugSnakemakeGenerator(
            benchmark_name="test",
            benchmark_version="1.0",
            benchmark_author="author",
        )

        output = StringIO()
        generator._write_gather_shell(output, gather_node)
        content = output.getvalue()

        # Should show execution preview
        assert 'echo "Would execute:"' in content
        assert (
            'echo "  cd {params.module_dir} && python3 {params.entrypoint} \\"'
            in content
        )
        assert f'echo "    --name {gather_node.module_id} \\"' in content
        assert 'echo "    --scores [files...] \\"' in content


@pytest.mark.short
class TestShellLineWriting:
    """Test the _write_shell_lines helper method."""

    def test_shell_lines_properly_formatted(self, basic_node):
        """Test that shell lines are properly formatted in shell: block."""
        generator = DebugSnakemakeGenerator(
            benchmark_name="test",
            benchmark_version="1.0",
            benchmark_author="author",
        )

        output = StringIO()
        lines = ['echo "test1"', 'echo "test2"', "touch output.txt"]
        generator._write_shell_lines(output, lines)
        content = output.getvalue()

        # Should have shell: directive
        assert "shell:" in content

        # All lines should be present
        assert 'echo "test1"' in content
        assert 'echo "test2"' in content
        assert "touch output.txt" in content


@pytest.mark.short
class TestIntegrationWithNodes:
    """Integration tests with various node configurations."""

    def test_mixed_node_types(self, basic_node, gather_node, node_with_parameters):
        """Test generator with mix of regular and gather nodes."""
        generator = DebugSnakemakeGenerator(
            benchmark_name="test",
            benchmark_version="1.0",
            benchmark_author="author",
        )

        # Should handle all node types without errors
        assert generator is not None

        # Test generating shell for each
        output1 = StringIO()
        generator._write_shell(output1, basic_node)
        assert len(output1.getvalue()) > 0

        output2 = StringIO()
        generator._write_gather_shell(output2, gather_node)
        assert len(output2.getvalue()) > 0

        output3 = StringIO()
        generator._write_shell(output3, node_with_parameters)
        assert len(output3.getvalue()) > 0

    def test_node_with_complex_input_mapping(self, basic_resolved_module):
        """Test node where input keys differ from original names."""
        node = ResolvedNode(
            id="complex_mapping",
            stage_id="stage",
            module_id="module",
            param_id="default",
            module=basic_resolved_module,
            inputs={
                "sanitized_input_0": "data/file1.csv",
                "sanitized_input_1": "data/file2.csv",
            },
            outputs=["output.txt"],
            input_name_mapping={
                "sanitized_input_0": "raw-data",
                "sanitized_input_1": "processed-data",
            },
        )

        generator = DebugSnakemakeGenerator(
            benchmark_name="test",
            benchmark_version="1.0",
            benchmark_author="author",
        )

        output = StringIO()
        generator._write_shell(output, node)
        content = output.getvalue()

        # Should map back to original names in execution preview
        assert 'echo "    --raw-data {input.sanitized_input_0} \\"' in content
        assert 'echo "    --processed-data {input.sanitized_input_1} \\"' in content

    def test_all_execution_paths_covered(
        self, basic_node, node_with_inputs, node_with_parameters, gather_node
    ):
        """Test that all code paths are exercised."""
        generator = DebugSnakemakeGenerator(
            benchmark_name="test",
            benchmark_version="1.0",
            benchmark_author="author",
        )

        # Write shell for regular nodes
        output1 = StringIO()
        generator._write_shell(output1, basic_node)
        assert len(output1.getvalue()) > 0

        output2 = StringIO()
        generator._write_shell(output2, node_with_inputs)
        assert len(output2.getvalue()) > 0

        output3 = StringIO()
        generator._write_shell(output3, node_with_parameters)
        assert len(output3.getvalue()) > 0

        # Write shell for gather node
        output4 = StringIO()
        generator._write_gather_shell(output4, gather_node)
        assert len(output4.getvalue()) > 0

        # All should produce unique content
        outputs = [
            output1.getvalue(),
            output2.getvalue(),
            output3.getvalue(),
            output4.getvalue(),
        ]
        assert len(set(outputs)) == 4  # All different
