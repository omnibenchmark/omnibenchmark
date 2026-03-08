"""
Edge case and corner case unit tests for SnakemakeGenerator.

Tests focus on boundary conditions, error paths, and unusual input combinations
that complement the main test_snakemake_gen.py coverage.

All tests marked as 'short' since they only test pure logic without I/O.
"""

import pytest
from pathlib import Path

from omnibenchmark.backend._metadata import save_metadata
from omnibenchmark.backend.snakemake import (
    SnakemakeGenerator,
    _bash_var,
    _make_human_name,
)
from omnibenchmark.model.params import Params
from omnibenchmark.model.benchmark import Resources, SoftwareBackendEnum, APIVersion
from omnibenchmark.model.resolved import (
    ResolvedNode,
    ResolvedModule,
    ResolvedEnvironment,
)


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
class TestBashVar:
    """Test _bash_var helper."""

    def test_dot_replaced(self):
        assert _bash_var("data.raw") == "_data_raw"

    def test_hyphen_replaced(self):
        assert _bash_var("my-input") == "_my_input"

    def test_plain_name_unchanged(self):
        assert _bash_var("counts") == "_counts"

    def test_multiple_dots(self):
        assert _bash_var("a.b.c") == "_a_b_c"

    def test_mixed(self):
        assert _bash_var("data.raw-v2") == "_data_raw_v2"


@pytest.mark.short
class TestMakeHumanNameEdgeCases:
    """Test edge cases in _make_human_name method."""

    def test_empty_dict(self, generator):
        """Test empty dict input."""
        result = _make_human_name({})
        assert result == ""

    def test_simple_params_dict(self, generator):
        """Test simple params dict."""
        params = {"alpha": 0.1, "beta": 0.5}
        result = _make_human_name(params)
        # Should create key-value pairs
        assert "alpha" in result
        assert "beta" in result

    def test_unsafe_characters_replaced(self, generator):
        """Test that unsafe characters are replaced."""
        params = {"key": "value/with\\unsafe:chars"}
        result = _make_human_name(params)
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
            output_path=output_path,
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
            output_path=output_path,
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
            output_path=output_path,
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
            output_path=output_path,
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
            output_path=output_path,
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
            output_path=output_path,
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
            output_path=output_path,
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
            output_path=output_path,
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
            output_path=output_path,
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
            output_path=output_path,
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
            output_path=output_path,
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

        save_metadata(benchmark_yaml, tmp_path, [])

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

        save_metadata(benchmark_yaml, tmp_path, nodes)

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
        save_metadata(benchmark_yaml, tmp_path, [])

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
            output_path=output_path,
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
            output_path=output_path,
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
            output_path=output_path,
        )

        content = output_path.read_text()
        # Should list all outputs
        assert "rule test_node:" in content


# ---------------------------------------------------------------------------
# Test _make_human_name long-name truncation (line 58)
# ---------------------------------------------------------------------------


@pytest.mark.short
class TestMakeHumanNameTruncation:
    """Test that names exceeding 255 chars are truncated with a hash suffix."""

    def test_long_name_truncated_with_hash(self):
        """Name > 255 chars falls back to [:246] + '_' + hash_short()."""
        # key=128 chars, value=128 chars -> "k-v" = 257 chars > 255
        long_key = "k" * 128
        long_val = "v" * 128
        params = Params({long_key: long_val})
        result = _make_human_name(params)
        assert len(result) <= 255
        # The 8-char short hash must appear at the end after the underscore
        assert result.endswith(params.hash_short())


# ---------------------------------------------------------------------------
# Test gather / collector node comments and shell dispatch (lines 112, 114, 167, 346-393)
# ---------------------------------------------------------------------------


@pytest.mark.short
class TestGatherCollectorNodes:
    """Test is_gather / is_collector branches in _write_node_rule."""

    def test_collector_node_gets_collector_comment(
        self, tmp_path, generator, minimal_module
    ):
        node = ResolvedNode(
            id="collector-node",
            stage_id="metrics",
            module_id="collector",
            param_id="default",
            module=minimal_module,
            is_collector=True,
            outputs=["metrics/out/.abc/result.json"],
        )
        output_path = tmp_path / "Snakefile"
        generator.generate_snakefile(nodes=[node], output_path=output_path)
        content = output_path.read_text()
        assert "# Type: Metric Collector" in content

    def test_gather_node_gets_gather_comment(self, tmp_path, generator, minimal_module):
        node = ResolvedNode(
            id="gather-node",
            stage_id="gather",
            module_id="gather",
            param_id="default",
            module=minimal_module,
            is_gather=True,
            outputs=["gather/out/.abc/result.json"],
        )
        output_path = tmp_path / "Snakefile"
        generator.generate_snakefile(nodes=[node], output_path=output_path)
        content = output_path.read_text()
        assert "# Type: Gather Stage" in content

    def test_gather_node_shell_uses_params_output_dir(
        self, tmp_path, generator, minimal_module
    ):
        """_write_gather_shell should use {params.output_dir} not $(dirname ...)."""
        node = ResolvedNode(
            id="gather-node",
            stage_id="gather",
            module_id="gather",
            param_id="default",
            module=minimal_module,
            is_gather=True,
            outputs=["gather/out/.abc/result.json"],
        )
        output_path = tmp_path / "Snakefile"
        generator.generate_snakefile(nodes=[node], output_path=output_path)
        content = output_path.read_text()
        assert "output_dir=" in content
        assert "{params.output_dir}" in content
        assert "$(dirname {output[0]})" not in content

    def test_gather_node_no_benchmark_directive(
        self, tmp_path, generator, minimal_module
    ):
        """Gather nodes must not emit benchmark: directive."""
        node = ResolvedNode(
            id="gather-node",
            stage_id="gather",
            module_id="gather",
            param_id="default",
            module=minimal_module,
            is_gather=True,
            outputs=["gather/out/.abc/result.json"],
        )
        output_path = tmp_path / "Snakefile"
        generator.generate_snakefile(nodes=[node], output_path=output_path)
        content = output_path.read_text()
        assert "benchmark:" not in content

    def test_gather_node_with_inputs_emits_bash_arrays(
        self, tmp_path, generator, minimal_module
    ):
        """Gather shell with inputs should build bash arrays per input name."""
        node = ResolvedNode(
            id="gather-node",
            stage_id="gather",
            module_id="gather",
            param_id="default",
            module=minimal_module,
            is_gather=True,
            inputs={"in_0": "data/D1/out.json", "in_1": "data/D2/out.json"},
            input_name_mapping={"in_0": "counts", "in_1": "counts"},
            outputs=["gather/out/.abc/result.json"],
        )
        output_path = tmp_path / "Snakefile"
        generator.generate_snakefile(nodes=[node], output_path=output_path)
        content = output_path.read_text()
        # bash array variable and loop must appear
        assert "_counts=()" in content
        assert "--counts" in content

    def test_gather_node_with_shebang_uses_direct_invocation(self, tmp_path, generator):
        """Gather shell with has_shebang=True must use $MODULE_DIR/entrypoint directly."""
        module = ResolvedModule(
            repository_url="https://github.com/test/repo.git",
            commit="abc123",
            module_dir=Path(".modules/abc123"),
            entrypoint=Path("run.py"),
            software_environment_id="test",
            has_shebang=True,
        )
        node = ResolvedNode(
            id="gather-node",
            stage_id="gather",
            module_id="gather",
            param_id="default",
            module=module,
            is_gather=True,
            outputs=["gather/out/.abc/result.json"],
        )
        output_path = tmp_path / "Snakefile"
        generator.generate_snakefile(nodes=[node], output_path=output_path)
        content = output_path.read_text()
        assert "$MODULE_DIR/{params.entrypoint}" in content

    def test_gather_node_with_parameters_emits_cli_args(
        self, tmp_path, generator, minimal_module
    ):
        """Gather shell with parameters should append {params.cli_args}."""
        params = Params({"k": "v"})
        node = ResolvedNode(
            id="gather-node",
            stage_id="gather",
            module_id="gather",
            param_id="p1",
            module=minimal_module,
            is_gather=True,
            parameters=params,
            outputs=["gather/out/.abc/result.json"],
        )
        output_path = tmp_path / "Snakefile"
        generator.generate_snakefile(nodes=[node], output_path=output_path)
        content = output_path.read_text()
        assert "{params.cli_args}" in content


# ---------------------------------------------------------------------------
# Test parameters branch (lines 136-137, 282-285, 332)
# ---------------------------------------------------------------------------


@pytest.mark.short
class TestNodeWithParameters:
    """Test that nodes with parameters emit cli_args and parameters.json lines."""

    def test_cli_args_emitted_in_params_section(
        self, tmp_path, generator, minimal_module
    ):
        params = Params({"alpha": 0.5, "beta": 1})
        node = ResolvedNode(
            id="test-node",
            stage_id="test",
            module_id="test",
            param_id="p1",
            module=minimal_module,
            parameters=params,
            outputs=["data/test/.abc/output.json"],
        )
        output_path = tmp_path / "Snakefile"
        generator.generate_snakefile(nodes=[node], output_path=output_path)
        content = output_path.read_text()
        assert "cli_args=" in content
        assert "--alpha" in content

    def test_params_json_written_in_shell(self, tmp_path, generator, minimal_module):
        params = Params({"k": "v"})
        node = ResolvedNode(
            id="test-node",
            stage_id="test",
            module_id="test",
            param_id="p1",
            module=minimal_module,
            parameters=params,
            outputs=["data/test/.abc/output.json"],
        )
        output_path = tmp_path / "Snakefile"
        generator.generate_snakefile(nodes=[node], output_path=output_path)
        content = output_path.read_text()
        assert "parameters.json" in content

    def test_cli_args_appended_to_command(self, tmp_path, generator, minimal_module):
        params = Params({"k": "v"})
        node = ResolvedNode(
            id="test-node",
            stage_id="test",
            module_id="test",
            param_id="p1",
            module=minimal_module,
            parameters=params,
            outputs=["data/test/.abc/output.json"],
        )
        output_path = tmp_path / "Snakefile"
        generator.generate_snakefile(nodes=[node], output_path=output_path)
        content = output_path.read_text()
        assert "{params.cli_args}" in content


# ---------------------------------------------------------------------------
# Test V0_5_0 benchmark path (line 148)
# ---------------------------------------------------------------------------


@pytest.mark.short
class TestBenchmarkDirective:
    """Test benchmark: directive paths for different API versions."""

    def test_v0_5_benchmark_uses_performance_txt(self, tmp_path, minimal_module):
        gen = SnakemakeGenerator(
            benchmark_name="Test",
            benchmark_version="1.0",
            benchmark_author="Author",
            api_version=APIVersion.V0_5_0,
        )
        node = ResolvedNode(
            id="test-node",
            stage_id="test",
            module_id="test",
            param_id="default",
            module=minimal_module,
            outputs=["data/test/.abc/output.json"],
        )
        output_path = tmp_path / "Snakefile"
        gen.generate_snakefile(nodes=[node], output_path=output_path)
        content = output_path.read_text()
        assert "performance.txt" in content

    def test_v0_4_benchmark_uses_dataset_name(self, tmp_path, minimal_module):
        gen = SnakemakeGenerator(
            benchmark_name="Test",
            benchmark_version="1.0",
            benchmark_author="Author",
            api_version=APIVersion.V0_4_0,
        )
        node = ResolvedNode(
            id="test-node",
            stage_id="test",
            module_id="test",
            param_id="default",
            module=minimal_module,
            outputs=["data/D1/.abc/D1_output.json"],
        )
        output_path = tmp_path / "Snakefile"
        gen.generate_snakefile(nodes=[node], output_path=output_path)
        content = output_path.read_text()
        assert "D1_performance.txt" in content


# ---------------------------------------------------------------------------
# Test environment directives (lines 196-199)
# ---------------------------------------------------------------------------


@pytest.mark.short
class TestEnvironmentDirectives:
    """Test that each backend type emits the correct Snakemake directive."""

    def _make_module(self, backend: SoftwareBackendEnum) -> ResolvedModule:
        env = ResolvedEnvironment(backend_type=backend, reference="envs/env.yaml")
        return ResolvedModule(
            repository_url="https://github.com/test/repo.git",
            commit="abc123",
            module_dir=Path(".modules/abc123"),
            entrypoint=Path("run.py"),
            software_environment_id="test",
            resolved_environment=env,
        )

    def test_apptainer_emits_container(self, tmp_path, generator):
        module = self._make_module(SoftwareBackendEnum.apptainer)
        node = ResolvedNode(
            id="test-node",
            stage_id="test",
            module_id="test",
            param_id="default",
            module=module,
            outputs=["output.json"],
        )
        output_path = tmp_path / "Snakefile"
        generator.generate_snakefile(nodes=[node], output_path=output_path)
        assert "container:" in output_path.read_text()

    def test_docker_emits_container(self, tmp_path, generator):
        module = self._make_module(SoftwareBackendEnum.docker)
        node = ResolvedNode(
            id="test-node",
            stage_id="test",
            module_id="test",
            param_id="default",
            module=module,
            outputs=["output.json"],
        )
        output_path = tmp_path / "Snakefile"
        generator.generate_snakefile(nodes=[node], output_path=output_path)
        assert "container:" in output_path.read_text()

    def test_envmodules_emits_envmodules(self, tmp_path, generator):
        module = self._make_module(SoftwareBackendEnum.envmodules)
        node = ResolvedNode(
            id="test-node",
            stage_id="test",
            module_id="test",
            param_id="default",
            module=module,
            outputs=["output.json"],
        )
        output_path = tmp_path / "Snakefile"
        generator.generate_snakefile(nodes=[node], output_path=output_path)
        assert "envmodules:" in output_path.read_text()


# ---------------------------------------------------------------------------
# Test resources directive with optional fields (lines 215-222)
# ---------------------------------------------------------------------------


@pytest.mark.short
class TestResourceDirectiveFields:
    """Test that optional resource fields are emitted when set."""

    def _node_with_resources(self, minimal_module, **kwargs) -> ResolvedNode:
        return ResolvedNode(
            id="test-node",
            stage_id="test",
            module_id="test",
            param_id="default",
            module=minimal_module,
            outputs=["output.json"],
            resources=Resources(**kwargs),
        )

    def test_mem_mb_emitted(self, tmp_path, generator, minimal_module):
        node = self._node_with_resources(minimal_module, cores=2, mem_mb=4096)
        output_path = tmp_path / "Snakefile"
        generator.generate_snakefile(nodes=[node], output_path=output_path)
        assert "mem_mb=4096" in output_path.read_text()

    def test_disk_mb_emitted(self, tmp_path, generator, minimal_module):
        node = self._node_with_resources(minimal_module, cores=2, disk_mb=8192)
        output_path = tmp_path / "Snakefile"
        generator.generate_snakefile(nodes=[node], output_path=output_path)
        assert "disk_mb=8192" in output_path.read_text()

    def test_runtime_emitted(self, tmp_path, generator, minimal_module):
        node = self._node_with_resources(minimal_module, cores=2, runtime=60)
        output_path = tmp_path / "Snakefile"
        generator.generate_snakefile(nodes=[node], output_path=output_path)
        assert "runtime=60" in output_path.read_text()

    def test_gpu_emitted(self, tmp_path, generator, minimal_module):
        node = self._node_with_resources(minimal_module, cores=2, gpu=1)
        output_path = tmp_path / "Snakefile"
        generator.generate_snakefile(nodes=[node], output_path=output_path)
        assert "nvidia_gpu=1" in output_path.read_text()


# ---------------------------------------------------------------------------
# Test shebang module (line 297)
# ---------------------------------------------------------------------------


@pytest.mark.short
class TestShebangedModule:
    """Test that has_shebang=True uses direct execution instead of interpreter."""

    def test_shebang_module_uses_direct_invocation(self, tmp_path, generator):
        module = ResolvedModule(
            repository_url="https://github.com/test/repo.git",
            commit="abc123",
            module_dir=Path(".modules/abc123"),
            entrypoint=Path("run.py"),
            software_environment_id="test",
            has_shebang=True,
        )
        node = ResolvedNode(
            id="test-node",
            stage_id="test",
            module_id="test",
            param_id="default",
            module=module,
            outputs=["data/test/.abc/output.json"],
        )
        output_path = tmp_path / "Snakefile"
        generator.generate_snakefile(nodes=[node], output_path=output_path)
        content = output_path.read_text()
        assert "./{params.entrypoint}" in content
        assert "python3" not in content


# ---------------------------------------------------------------------------
# Test output_dir param in generated Snakefile (new behaviour)
# ---------------------------------------------------------------------------


@pytest.mark.short
class TestOutputDirParam:
    """Ensure output_dir is emitted as an explicit param and used in shell."""

    def test_output_dir_in_params_section(self, tmp_path, generator, minimal_module):
        node = ResolvedNode(
            id="datasets-pbmc3k-default",
            stage_id="datasets",
            module_id="pbmc3k",
            param_id="default",
            module=minimal_module,
            outputs=["data/pbmc3k/.d41d8cd9/pbmc3k.h5ad"],
        )
        output_path = tmp_path / "Snakefile"
        generator.generate_snakefile(nodes=[node], output_path=output_path)
        content = output_path.read_text()
        assert 'output_dir="data/pbmc3k/.d41d8cd9"' in content

    def test_shell_uses_params_output_dir_not_dirname(
        self, tmp_path, generator, minimal_module
    ):
        node = ResolvedNode(
            id="test-node",
            stage_id="test",
            module_id="test",
            param_id="default",
            module=minimal_module,
            outputs=["data/test/.abc/output.json"],
        )
        output_path = tmp_path / "Snakefile"
        generator.generate_snakefile(nodes=[node], output_path=output_path)
        content = output_path.read_text()
        assert "{params.output_dir}" in content
        assert "$(dirname {output[0]})" not in content

    def test_node_with_no_outputs_has_no_output_dir(
        self, tmp_path, generator, minimal_module
    ):
        node = ResolvedNode(
            id="test-node",
            stage_id="test",
            module_id="test",
            param_id="default",
            module=minimal_module,
            outputs=[],
        )
        output_path = tmp_path / "Snakefile"
        generator.generate_snakefile(nodes=[node], output_path=output_path)
        content = output_path.read_text()
        assert "output_dir=" not in content
