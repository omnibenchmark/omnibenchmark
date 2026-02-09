"""Tests for Snakemake generation of gather stage nodes."""

import io
from pathlib import Path


from omnibenchmark.backend.snakemake_gen import SnakemakeGenerator
from omnibenchmark.model.resolved import ResolvedModule, ResolvedNode, TemplateContext


def _make_module(**overrides):
    """Create a minimal ResolvedModule for testing."""
    defaults = dict(
        repository_url="https://github.com/test/repo",
        commit="a" * 40,
        module_dir=Path(".modules/repo/aaa"),
        entrypoint=Path("run.py"),
        software_environment_id="default",
        resolved_environment=None,
        has_shebang=False,
        interpreter="python3",
        dirty=False,
    )
    defaults.update(overrides)
    return ResolvedModule(**defaults)


def _make_gather_node(**overrides):
    """Create a minimal gather ResolvedNode for testing."""
    module = overrides.pop("module", _make_module())
    defaults = dict(
        id="gather_aggregate-collector_mod.default",
        stage_id="aggregate",
        module_id="collector_mod",
        param_id=".default",
        module=module,
        parameters=None,
        inputs={
            "input_0": "methods/mod_a/.default/result.csv",
            "input_1": "methods/mod_b/.default/result.csv",
        },
        outputs=["aggregate/collector_mod/.default/summary.csv"],
        template_context=TemplateContext(
            provides={},
            module_attrs={"id": "collector_mod", "stage": "aggregate"},
        ),
        input_name_mapping={
            "input_0": "method",
            "input_1": "method",
        },
    )
    defaults.update(overrides)
    return ResolvedNode(**defaults)


def _make_regular_node(**overrides):
    """Create a minimal regular (map) ResolvedNode for testing."""
    module = overrides.pop("module", _make_module())
    defaults = dict(
        id="methods-mod_a-.default",
        stage_id="methods",
        module_id="mod_a",
        param_id=".default",
        module=module,
        parameters=None,
        inputs={"data_matrix": "data/dataset1/.default/matrix.csv"},
        outputs=["methods/mod_a/.default/result.csv"],
        template_context=TemplateContext(
            provides={"dataset": "dataset1"},
            module_attrs={"id": "mod_a", "stage": "methods"},
        ),
        input_name_mapping={"data_matrix": "data.matrix"},
    )
    defaults.update(overrides)
    return ResolvedNode(**defaults)


def _make_collector_node(**overrides):
    """Create a minimal collector ResolvedNode for testing."""
    module = overrides.pop("module", _make_module())
    defaults = dict(
        id="_collector_metrics-metric_mod-.default",
        stage_id="_collector_metrics",
        module_id="metric_mod",
        param_id=".default",
        module=module,
        parameters=None,
        inputs={
            "input_0": "methods/mod_a/.default/result.csv",
            "input_1": "methods/mod_b/.default/result.csv",
        },
        outputs=["_collector_metrics/metric_mod/.default/metrics.json"],
        template_context=TemplateContext(
            provides={},
            module_attrs={"id": "metric_mod", "stage": "_collector_metrics"},
        ),
        input_name_mapping={
            "input_0": "result",
            "input_1": "result",
        },
    )
    defaults.update(overrides)
    return ResolvedNode(**defaults)


def _generate_rule(node, debug_mode=True):
    """Generate the Snakemake rule text for a single node."""
    gen = SnakemakeGenerator("test_bench", "1.0", "tester")
    buf = io.StringIO()
    gen._write_node_rule(buf, node, debug_mode=debug_mode)
    return buf.getvalue()


class TestGatherNodeRuleGeneration:
    """Tests for gather node Snakemake rule generation."""

    def test_gather_node_has_type_comment(self):
        """Gather nodes should have '# Type: Gather Stage' comment."""
        node = _make_gather_node()
        rule = _generate_rule(node)
        assert "# Type: Gather Stage" in rule

    def test_regular_node_has_no_type_comment(self):
        """Regular nodes should not have type comments."""
        node = _make_regular_node()
        rule = _generate_rule(node)
        assert "# Type:" not in rule

    def test_collector_node_has_collector_type(self):
        """Collector nodes should have '# Type: Metric Collector' comment."""
        node = _make_collector_node()
        rule = _generate_rule(node)
        assert "# Type: Metric Collector" in rule

    def test_gather_node_no_benchmark_directive(self):
        """Gather nodes should NOT have a benchmark: directive."""
        node = _make_gather_node()
        rule = _generate_rule(node)
        assert "    benchmark:" not in rule
        assert "clustbench_performance" not in rule

    def test_regular_node_has_benchmark_directive(self):
        """Regular nodes SHOULD have a benchmark: directive."""
        node = _make_regular_node()
        rule = _generate_rule(node)
        assert "    benchmark:" in rule
        assert "clustbench_performance" in rule

    def test_collector_node_no_benchmark_directive(self):
        """Collector nodes should NOT have a benchmark: directive."""
        node = _make_collector_node()
        rule = _generate_rule(node)
        assert "    benchmark:" not in rule

    def test_gather_node_has_named_inputs(self):
        """Gather node inputs should be written as named entries."""
        node = _make_gather_node()
        rule = _generate_rule(node)
        assert 'input_0="methods/mod_a/.default/result.csv"' in rule
        assert 'input_1="methods/mod_b/.default/result.csv"' in rule

    def test_gather_node_has_outputs(self):
        """Gather node outputs should be written correctly."""
        node = _make_gather_node()
        rule = _generate_rule(node)
        assert '"aggregate/collector_mod/.default/summary.csv"' in rule

    def test_gather_node_has_log_directive(self):
        """Gather nodes should have a log: directive."""
        node = _make_gather_node()
        rule = _generate_rule(node)
        assert "    log:" in rule
        assert ".logs/" in rule

    def test_gather_node_has_resources(self):
        """Gather nodes should have resources directive."""
        node = _make_gather_node()
        rule = _generate_rule(node)
        assert "    resources:" in rule


class TestGatherDebugShell:
    """Tests for gather node debug shell output."""

    def test_debug_shell_labels_as_gather(self):
        """Debug shell should echo 'GATHER STAGE' not 'METRIC COLLECTOR'."""
        node = _make_gather_node()
        rule = _generate_rule(node, debug_mode=True)
        assert "GATHER STAGE" in rule
        assert "METRIC COLLECTOR" not in rule

    def test_collector_debug_shell_labels_as_collector(self):
        """Collector debug shell should still echo 'METRIC COLLECTOR'."""
        node = _make_collector_node()
        rule = _generate_rule(node, debug_mode=True)
        assert "METRIC COLLECTOR" in rule
        assert "GATHER STAGE" not in rule

    def test_debug_shell_shows_module_id(self):
        """Debug shell should include the module ID."""
        node = _make_gather_node()
        rule = _generate_rule(node, debug_mode=True)
        assert "collector_mod" in rule

    def test_debug_shell_creates_output_dirs(self):
        """Debug shell should create output directories."""
        node = _make_gather_node()
        rule = _generate_rule(node, debug_mode=True)
        assert "mkdir -p" in rule

    def test_debug_shell_touches_outputs(self):
        """Debug shell should touch output files."""
        node = _make_gather_node()
        rule = _generate_rule(node, debug_mode=True)
        assert "touch" in rule


class TestGatherExecShell:
    """Tests for gather node execution shell output."""

    def test_exec_shell_creates_output_dirs(self):
        """Exec shell should create output directories."""
        node = _make_gather_node()
        rule = _generate_rule(node, debug_mode=False)
        assert "mkdir -p" in rule

    def test_exec_shell_uses_module_dir(self):
        """Exec shell should reference the module directory."""
        node = _make_gather_node()
        rule = _generate_rule(node, debug_mode=False)
        assert "MODULE_DIR" in rule

    def test_exec_shell_passes_output_dir(self):
        """Exec shell should pass --output_dir."""
        node = _make_gather_node()
        rule = _generate_rule(node, debug_mode=False)
        assert "--output_dir" in rule

    def test_exec_shell_passes_name(self):
        """Exec shell should pass --name with module_id."""
        node = _make_gather_node()
        rule = _generate_rule(node, debug_mode=False)
        assert "--name collector_mod" in rule

    def test_exec_shell_passes_inputs_by_name(self):
        """Exec shell should pass input files grouped by their original name."""
        node = _make_gather_node()
        rule = _generate_rule(node, debug_mode=False)
        assert "--method" in rule


class TestGatherInFullSnakefile:
    """Tests for gather nodes in a complete generated Snakefile."""

    def test_full_snakefile_includes_gather_and_regular(self):
        """A Snakefile with both regular and gather nodes should include both."""
        regular = _make_regular_node()
        gather = _make_gather_node()
        gen = SnakemakeGenerator("test_bench", "1.0", "tester")

        _buf = io.StringIO()
        with open("/dev/null", "w"):
            pass

        # Generate to a temp file
        import tempfile

        with tempfile.NamedTemporaryFile(mode="w", suffix=".smk", delete=False) as f:
            tmp_path = Path(f.name)

        try:
            gen.generate_snakefile(
                nodes=[regular, gather],
                collectors=[],
                output_path=tmp_path,
                debug_mode=True,
            )
            content = tmp_path.read_text()

            # Both rules present
            assert "rule methods_mod_a__default" in content
            assert "rule gather_aggregate_collector_mod_default" in content

            # All rule should include outputs from both
            assert "rule all:" in content
            assert "methods/mod_a/.default/result.csv" in content
            assert "aggregate/collector_mod/.default/summary.csv" in content
        finally:
            tmp_path.unlink()

    def test_gather_node_outputs_in_all_rule(self):
        """Gather node outputs should appear in the 'all' rule."""
        gather = _make_gather_node()
        gen = SnakemakeGenerator("test_bench", "1.0", "tester")

        _buf = io.StringIO()
        gen._write_all_rule(_buf, [gather], [])
        content = _buf.getvalue()

        assert "aggregate/collector_mod/.default/summary.csv" in content

    def test_sanitized_rule_name(self):
        """Gather node IDs should be sanitized to valid Snakemake rule names."""
        gen = SnakemakeGenerator("test_bench", "1.0", "tester")
        name = gen._sanitize_rule_name("gather_aggregate-collector_mod.default")
        # Hyphens and dots become underscores
        assert "-" not in name
        assert "." not in name
        assert name == "gather_aggregate_collector_mod_default"
