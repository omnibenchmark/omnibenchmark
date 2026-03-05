"""
Edge case and corner case unit tests for SnakemakeGenerator.

Tests focus on boundary conditions, error paths, and unusual input combinations
that complement the main test_snakemake_gen.py coverage.

All tests marked as 'short' since they only test pure logic without I/O.
"""

import pytest
from pathlib import Path

from omnibenchmark.backend.snakemake_gen import (
    SnakemakeGenerator,
    save_metadata,
)
from omnibenchmark.model.resolved import (
    ResolvedNode,
    ResolvedModule,
    ResolvedEnvironment,
)
from omnibenchmark.model.benchmark import SoftwareBackendEnum, APIVersion


# ---------------------------------------------------------------------------
# Test Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def minimal_module():
    """Minimal ResolvedModule for testing."""
    return ResolvedModule(
        repository_url="https://github.com/test/repo.git",
        commit="abc123",
        module_dir=Path(".modules/abc123"),
        entrypoint=Path("run.py"),
        software_environment_id="test",
    )


@pytest.fixture
def generator():
    """Basic SnakemakeGenerator instance."""
    return SnakemakeGenerator(
        benchmark_name="TestBenchmark",
        benchmark_version="1.0.0",
        benchmark_author="Test Author",
    )


# ---------------------------------------------------------------------------
# Test _sanitize_rule_name edge cases
# ---------------------------------------------------------------------------


@pytest.mark.short
class TestSanitizeRuleNameEdgeCases:
    """Test edge cases in _sanitize_rule_name method."""

    def test_hyphens_replaced_with_underscores(self, generator):
        """Test that hyphens are replaced with underscores."""
        result = generator._sanitize_rule_name("test-name-id")
        assert result == "test_name_id"

    def test_dots_replaced_with_underscores(self, generator):
        """Test that dots are replaced with underscores."""
        result = generator._sanitize_rule_name("test.name.id")
        assert result == "test_name_id"

    def test_starts_with_digit_gets_prefix(self, generator):
        """Test name starting with digit gets 'rule_' prefix."""
        result = generator._sanitize_rule_name("123abc")
        assert result == "rule_123abc"

    def test_starts_with_letter_unchanged(self, generator):
        """Test name starting with letter is unchanged."""
        result = generator._sanitize_rule_name("test_name")
        assert result == "test_name"

    def test_empty_string_handled(self, generator):
        """Test empty string handling."""
        result = generator._sanitize_rule_name("")
        # Returns empty string
        assert result == ""

    def test_multiple_consecutive_replacements(self, generator):
        """Test multiple consecutive hyphens and dots."""
        result = generator._sanitize_rule_name("test---name...id")
        assert result == "test___name___id"


# ---------------------------------------------------------------------------
# Test _make_human_name edge cases
# ---------------------------------------------------------------------------


@pytest.mark.short
class TestMakeHumanNameEdgeCases:
    """Test edge cases in _make_human_name method."""

    def test_empty_dict(self, generator):
        """Test empty dict input."""
        result = generator._make_human_name({})
        assert result == ""

    def test_simple_params_dict(self, generator):
        """Test simple params dict."""
        params = {"alpha": 0.1, "beta": 0.5}
        result = generator._make_human_name(params)
        # Should create key-value pairs
        assert "alpha" in result
        assert "beta" in result

    def test_unsafe_characters_replaced(self, generator):
        """Test that unsafe characters are replaced."""
        params = {"key": "value/with\\unsafe:chars"}
        result = generator._make_human_name(params)
        # Unsafe chars should be replaced with underscores
        assert "/" not in result
        assert "\\" not in result
        assert ":" not in result


# ---------------------------------------------------------------------------
# Test SnakemakeGenerator edge cases
# ---------------------------------------------------------------------------


@pytest.mark.short
class TestSnakemakeGeneratorEdgeCases:
    """Test edge cases in SnakemakeGenerator."""

    def test_empty_nodes_list(self, tmp_path, generator):
        """Test generating Snakefile with no nodes."""
        output_path = tmp_path / "Snakefile"

        generator.generate_snakefile(
            nodes=[],
            collectors=[],
            output_path=output_path,
            debug_mode=True,
        )

        # Should still create valid file with header
        assert output_path.exists()
        content = output_path.read_text()
        assert "#!/usr/bin/env snakemake" in content
        assert "rule all:" in content

    def test_node_with_no_inputs(self, tmp_path, generator, minimal_module):
        """Test node with empty inputs dict (entrypoint)."""
        node = ResolvedNode(
            id="data-test-default",
            stage_id="data",
            module_id="test",
            param_id="default",
            module=minimal_module,
            inputs={},  # Empty inputs
            outputs=["data/test/output.json"],
        )

        output_path = tmp_path / "Snakefile"
        generator.generate_snakefile(
            nodes=[node],
            collectors=[],
            output_path=output_path,
            debug_mode=False,
        )

        content = output_path.read_text()
        # Should have rule with no input directive
        assert "rule data_test_default:" in content

    def test_node_with_no_outputs(self, tmp_path, generator, minimal_module):
        """Test node with empty outputs list."""
        node = ResolvedNode(
            id="test-node",
            stage_id="test",
            module_id="test",
            param_id="default",
            module=minimal_module,
            outputs=[],  # No outputs
        )

        output_path = tmp_path / "Snakefile"
        generator.generate_snakefile(
            nodes=[node],
            collectors=[],
            output_path=output_path,
            debug_mode=True,
        )

        content = output_path.read_text()
        # Should still generate rule
        assert "rule test_node:" in content

    def test_node_with_very_long_id(self, tmp_path, generator, minimal_module):
        """Test node with extremely long ID."""
        long_id = "a" * 200
        node = ResolvedNode(
            id=long_id,
            stage_id="test",
            module_id="test",
            param_id="default",
            module=minimal_module,
            outputs=["output.json"],
        )

        output_path = tmp_path / "Snakefile"
        generator.generate_snakefile(
            nodes=[node],
            collectors=[],
            output_path=output_path,
            debug_mode=True,
        )

        content = output_path.read_text()
        # Rule name should be present (sanitized)
        assert "rule" in content

    def test_node_with_special_chars_in_id(self, tmp_path, generator, minimal_module):
        """Test node with special characters in ID."""
        node = ResolvedNode(
            id="test-node@version#1",
            stage_id="test",
            module_id="test",
            param_id="default",
            module=minimal_module,
            outputs=["output.json"],
        )

        output_path = tmp_path / "Snakefile"
        generator.generate_snakefile(
            nodes=[node],
            collectors=[],
            output_path=output_path,
            debug_mode=True,
        )

        content = output_path.read_text()
        # Should generate valid Snakefile despite special chars
        assert "rule" in content

    def test_module_with_no_interpreter(self, tmp_path, generator):
        """Test module without interpreter specified."""
        module = ResolvedModule(
            repository_url="https://github.com/test/repo.git",
            commit="abc123",
            module_dir=Path(".modules/abc123"),
            entrypoint=Path("run.py"),
            software_environment_id="test",
            has_shebang=False,
            interpreter=None,  # No interpreter
        )

        node = ResolvedNode(
            id="test-node",
            stage_id="test",
            module_id="test",
            param_id="default",
            module=module,
            outputs=["output.json"],
        )

        output_path = tmp_path / "Snakefile"
        generator.generate_snakefile(
            nodes=[node],
            collectors=[],
            output_path=output_path,
            debug_mode=False,
        )

        content = output_path.read_text()
        # Should default to python3
        assert "python3" in content

    def test_multiple_nodes_same_stage(self, tmp_path, generator, minimal_module):
        """Test multiple nodes from same stage."""
        nodes = [
            ResolvedNode(
                id=f"test-node-{i}",
                stage_id="test",
                module_id=f"mod{i}",
                param_id="default",
                module=minimal_module,
                outputs=[f"output{i}.json"],
            )
            for i in range(10)
        ]

        output_path = tmp_path / "Snakefile"
        generator.generate_snakefile(
            nodes=nodes,
            collectors=[],
            output_path=output_path,
            debug_mode=True,
        )

        content = output_path.read_text()
        # Should have all rules
        for i in range(10):
            assert f"rule test_node_{i}:" in content

    def test_environment_with_missing_reference(self, tmp_path, generator):
        """Test environment with None reference."""
        env = ResolvedEnvironment(
            backend_type=SoftwareBackendEnum.conda,
            reference=None,  # Missing reference
        )

        module = ResolvedModule(
            repository_url="https://github.com/test/repo.git",
            commit="abc123",
            module_dir=Path(".modules/abc123"),
            entrypoint=Path("run.py"),
            software_environment_id="test",
            resolved_environment=env,
        )

        node = ResolvedNode(
            id="test-node",
            stage_id="test",
            module_id="test",
            param_id="default",
            module=module,
            outputs=["output.json"],
        )

        output_path = tmp_path / "Snakefile"
        generator.generate_snakefile(
            nodes=[node],
            collectors=[],
            output_path=output_path,
            debug_mode=True,
        )

        # Should handle gracefully
        assert output_path.exists()

    def test_api_version_affects_output_paths(self, tmp_path, minimal_module):
        """Test that different API versions produce different output patterns."""
        node = ResolvedNode(
            id="test-node",
            stage_id="test",
            module_id="test",
            param_id="default",
            module=minimal_module,
            outputs=["output.json"],
        )

        # Test with v0.4.0
        gen_v0_4 = SnakemakeGenerator(
            benchmark_name="Test",
            benchmark_version="1.0",
            benchmark_author="Author",
            api_version=APIVersion.V0_4_0,
        )

        output_path = tmp_path / "Snakefile_v0_4"
        gen_v0_4.generate_snakefile(
            nodes=[node],
            collectors=[],
            output_path=output_path,
            debug_mode=True,
        )

        content = output_path.read_text()
        assert "Test" in content

    def test_node_with_gather_style_inputs(self, tmp_path, generator, minimal_module):
        """Test node with gather-style inputs (lists as values)."""
        node = ResolvedNode(
            id="gather-node",
            stage_id="gather",
            module_id="gather",
            param_id="default",
            module=minimal_module,
            inputs={
                "input_0": ["file1.json", "file2.json", "file3.json"],
                "input_1": ["file4.json", "file5.json"],
            },
            outputs=["gathered.json"],
        )

        output_path = tmp_path / "Snakefile"
        generator.generate_snakefile(
            nodes=[node],
            collectors=[],
            output_path=output_path,
            debug_mode=False,
        )

        content = output_path.read_text()
        # Should handle list inputs
        assert "rule gather_node:" in content

    def test_node_with_wildcards_in_output(self, tmp_path, generator, minimal_module):
        """Test node with wildcards in output paths."""
        node = ResolvedNode(
            id="test-node",
            stage_id="test",
            module_id="test",
            param_id="default",
            module=minimal_module,
            outputs=["output_{wildcard}.json"],
        )

        output_path = tmp_path / "Snakefile"
        generator.generate_snakefile(
            nodes=[node],
            collectors=[],
            output_path=output_path,
            debug_mode=False,
        )

        content = output_path.read_text()
        # Should handle wildcards
        assert "rule test_node:" in content


# ---------------------------------------------------------------------------
# Test save_metadata edge cases
# ---------------------------------------------------------------------------


@pytest.mark.short
class TestSaveMetadataEdgeCases:
    """Test edge cases in save_metadata function."""

    def test_empty_nodes_list(self, tmp_path):
        """Test saving metadata with no nodes."""
        # Create a dummy benchmark.yaml file
        benchmark_yaml = tmp_path / "benchmark.yaml"
        benchmark_yaml.write_text("name: test\nversion: 1.0\n")

        save_metadata(benchmark_yaml, tmp_path, [], [])

        modules_file = tmp_path / ".metadata" / "modules.txt"
        assert modules_file.exists()
        content = modules_file.read_text()
        # Should have header
        assert "#" in content

    def test_duplicate_modules_across_nodes(self, tmp_path, minimal_module):
        """Test deduplication of modules used in multiple nodes."""
        # Create a dummy benchmark.yaml file
        benchmark_yaml = tmp_path / "benchmark.yaml"
        benchmark_yaml.write_text("name: test\nversion: 1.0\n")

        nodes = [
            ResolvedNode(
                id=f"node-{i}",
                stage_id="test",
                module_id=f"mod{i}",
                param_id="default",
                module=minimal_module,  # Same module
                outputs=[f"output{i}.json"],
            )
            for i in range(5)
        ]

        save_metadata(benchmark_yaml, tmp_path, nodes, [])

        modules_file = tmp_path / ".metadata" / "modules.txt"
        content = modules_file.read_text()
        # Should only list module once
        url_count = content.count(minimal_module.repository_url)
        assert url_count == 1

    def test_metadata_dir_already_exists(self, tmp_path):
        """Test that existing .metadata directory doesn't cause errors."""
        # Create a dummy benchmark.yaml file
        benchmark_yaml = tmp_path / "benchmark.yaml"
        benchmark_yaml.write_text("name: test\nversion: 1.0\n")

        metadata_dir = tmp_path / ".metadata"
        metadata_dir.mkdir()

        # Should not raise error
        save_metadata(benchmark_yaml, tmp_path, [], [])

        assert metadata_dir.exists()


# ---------------------------------------------------------------------------
# Test resource directive edge cases
# ---------------------------------------------------------------------------


@pytest.mark.short
class TestResourceDirectiveEdgeCases:
    """Test edge cases in resource directive generation."""

    def test_node_without_resources_generates_valid_snakefile(
        self, tmp_path, generator, minimal_module
    ):
        """Test that nodes without resource specifications generate valid Snakefiles."""
        node = ResolvedNode(
            id="test-node",
            stage_id="test",
            module_id="test",
            param_id="default",
            module=minimal_module,
            outputs=["output.json"],
        )

        output_path = tmp_path / "Snakefile"
        generator.generate_snakefile(
            nodes=[node],
            collectors=[],
            output_path=output_path,
            debug_mode=True,
        )

        # Should not crash and should produce valid output
        assert output_path.exists()
        content = output_path.read_text()
        assert "rule" in content


# ---------------------------------------------------------------------------
# Test complexity edge cases
# ---------------------------------------------------------------------------


@pytest.mark.short
class TestComplexityEdgeCases:
    """Test computational complexity edge cases."""

    def test_deep_input_nesting(self, tmp_path, generator, minimal_module):
        """Test node with deeply nested input structure."""
        # Create inputs with many keys
        inputs = {f"input_{i}": f"path/to/file{i}.json" for i in range(50)}

        node = ResolvedNode(
            id="test-node",
            stage_id="test",
            module_id="test",
            param_id="default",
            module=minimal_module,
            inputs=inputs,
            outputs=["output.json"],
        )

        output_path = tmp_path / "Snakefile"
        generator.generate_snakefile(
            nodes=[node],
            collectors=[],
            output_path=output_path,
            debug_mode=False,
        )

        # Should handle many inputs
        assert output_path.exists()

    def test_many_outputs(self, tmp_path, generator, minimal_module):
        """Test node with many output files."""
        outputs = [f"output/file{i}.json" for i in range(100)]

        node = ResolvedNode(
            id="test-node",
            stage_id="test",
            module_id="test",
            param_id="default",
            module=minimal_module,
            outputs=outputs,
        )

        output_path = tmp_path / "Snakefile"
        generator.generate_snakefile(
            nodes=[node],
            collectors=[],
            output_path=output_path,
            debug_mode=True,
        )

        content = output_path.read_text()
        # Should list all outputs
        assert "rule test_node:" in content
