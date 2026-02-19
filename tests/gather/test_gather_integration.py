"""Integration tests for the full gather pipeline.

Tests the complete gather flow from YAML parsing through model validation,
ordering validation, node expansion, and Snakemake generation — without
requiring real git repos or module resolution.
"""

import io
import textwrap
from pathlib import Path

import pytest

from omnibenchmark.backend.snakemake_gen import SnakemakeGenerator
from omnibenchmark.cli.run import (
    GatherOrderingError,
    _build_gather_inputs,
    _validate_gather_ordering,
)
from omnibenchmark.model.benchmark import Benchmark
from omnibenchmark.model.resolved import ResolvedModule, ResolvedNode


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_resolved_module(**overrides):
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


FAKE_MODULE = _make_resolved_module()


# ---------------------------------------------------------------------------
# 1. Full pipeline: YAML → Model → Validation → Expansion → Snakemake
# ---------------------------------------------------------------------------


class TestGatherPipelineIntegration:
    """End-to-end integration: parse YAML, validate, expand nodes, generate Snakefile."""

    YAML_CONFIG = textwrap.dedent("""\
        id: GatherIntegrationBenchmark
        description: Integration test for provides/gather
        version: "1.0"
        benchmarker: "Test"
        benchmark_yaml_spec: "0.3"
        software_backend: host

        software_environments:
          host:
            description: "Host"

        stages:
          - id: data
            modules:
              - id: D1
                name: "Dataset 1"
                software_environment: "host"
                repository:
                  url: https://github.com/test/data
                  commit: aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
            outputs:
              - id: data.raw
                path: "{dataset}_data.csv"

          - id: method_a
            provides:
              method_output: method_a.result
            modules:
              - id: MA
                name: "Method A"
                software_environment: "host"
                repository:
                  url: https://github.com/test/method_a
                  commit: bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb
            inputs:
              - data.raw
            outputs:
              - id: method_a.result
                path: "{dataset}_result.csv"

          - id: method_b
            provides:
              method_output: method_b.result
            modules:
              - id: MB
                name: "Method B"
                software_environment: "host"
                repository:
                  url: https://github.com/test/method_b
                  commit: cccccccccccccccccccccccccccccccccccccccc
            inputs:
              - data.raw
            outputs:
              - id: method_b.result
                path: "{dataset}_result.csv"

          - id: aggregate
            modules:
              - id: AGG
                name: "Aggregator"
                software_environment: "host"
                repository:
                  url: https://github.com/test/aggregator
                  commit: dddddddddddddddddddddddddddddddddddddd
            inputs:
              - gather: method_output
            outputs:
              - id: aggregate.summary
                path: summary.csv
    """)

    def _parse_benchmark(self):
        """Parse the YAML and return a Benchmark model."""
        return Benchmark.from_yaml(self.YAML_CONFIG)

    def test_yaml_parses_successfully(self):
        """The YAML with provides/gather should parse into a valid Benchmark."""
        bm = self._parse_benchmark()
        assert len(bm.stages) == 4
        assert bm.stages[0].id == "data"
        assert bm.stages[3].id == "aggregate"

    def test_provides_parsed_correctly(self):
        """Stages with 'provides' should have the labels parsed."""
        bm = self._parse_benchmark()
        assert bm.stages[1].provides == {"method_output": "method_a.result"}
        assert bm.stages[2].provides == {"method_output": "method_b.result"}
        assert bm.stages[0].provides is None
        assert bm.stages[3].provides is None

    def test_gather_input_parsed_correctly(self):
        """The aggregate stage should have a GatherInput."""
        bm = self._parse_benchmark()
        agg = bm.stages[3]
        assert agg.is_gather_stage()
        assert agg.get_gather_labels() == ["method_output"]

    def test_provider_stages_found(self):
        """get_provider_stages should return method_a and method_b."""
        bm = self._parse_benchmark()
        providers = bm.get_provider_stages("method_output")
        provider_ids = [s.id for s in providers]
        assert "method_a" in provider_ids
        assert "method_b" in provider_ids
        assert len(provider_ids) == 2

    def test_ordering_validation_passes(self):
        """Providers before gatherer — ordering should pass."""
        bm = self._parse_benchmark()
        # Should not raise
        _validate_gather_ordering(bm.stages)

    def test_ordering_validation_fails_when_reversed(self):
        """Swap aggregate before providers — should fail."""
        bm = self._parse_benchmark()
        # Put aggregate first
        reordered = [bm.stages[3], bm.stages[0], bm.stages[1], bm.stages[2]]
        with pytest.raises(GatherOrderingError):
            _validate_gather_ordering(reordered)

    def test_node_expansion_builds_gather_inputs(self):
        """_build_gather_inputs should collect outputs from provider nodes."""
        bm = self._parse_benchmark()

        # Build output_to_nodes as the expansion loop would
        output_to_nodes = {
            "method_a.result": [
                ("method_a-MA-.default-D1", "method_a/MA/.default/D1_result.csv")
            ],
            "method_b.result": [
                ("method_b-MB-.default-D1", "method_b/MB/.default/D1_result.csv")
            ],
        }

        inputs, name_mapping = _build_gather_inputs(
            ["method_output"], [], bm, output_to_nodes
        )

        assert len(inputs) == 2
        assert inputs["input_0"] == "method_a/MA/.default/D1_result.csv"
        assert inputs["input_1"] == "method_b/MB/.default/D1_result.csv"
        assert name_mapping["input_0"] == "method_output"
        assert name_mapping["input_1"] == "method_output"

    def test_snakemake_generation_for_gather_node(self):
        """A gather ResolvedNode should produce a valid Snakemake rule."""
        gather_node = ResolvedNode(
            id="gather_aggregate-AGG-.default",
            stage_id="aggregate",
            module_id="AGG",
            param_id=".default",
            module=FAKE_MODULE,
            inputs={
                "input_0": "method_a/MA/.default/D1_result.csv",
                "input_1": "method_b/MB/.default/D1_result.csv",
            },
            outputs=["aggregate/AGG/.default/summary.csv"],
            input_name_mapping={
                "input_0": "method_output",
                "input_1": "method_output",
            },
        )

        gen = SnakemakeGenerator("GatherIntegrationBenchmark", "1.0", "Test")
        buf = io.StringIO()
        gen._write_node_rule(buf, gather_node, debug_mode=True)
        rule_text = buf.getvalue()

        # Verify structure
        assert "# Type: Gather Stage" in rule_text
        assert "rule gather_aggregate_AGG__default:" in rule_text
        assert 'input_0="method_a/MA/.default/D1_result.csv"' in rule_text
        assert 'input_1="method_b/MB/.default/D1_result.csv"' in rule_text
        assert '"aggregate/AGG/.default/summary.csv"' in rule_text
        assert "benchmark:" not in rule_text
        assert "GATHER STAGE" in rule_text

    def test_full_snakefile_with_gather(self):
        """Generate a complete Snakefile with data, method, and gather nodes."""
        data_node = ResolvedNode(
            id="data-D1-.default",
            stage_id="data",
            module_id="D1",
            param_id=".default",
            module=FAKE_MODULE,
            outputs=["data/D1/.default/D1_data.csv"],
        )
        method_a_node = ResolvedNode(
            id="method_a-MA-.default-D1",
            stage_id="method_a",
            module_id="MA",
            param_id=".default",
            module=FAKE_MODULE,
            inputs={"data_raw": "data/D1/.default/D1_data.csv"},
            outputs=["method_a/MA/.default/D1_result.csv"],
            input_name_mapping={"data_raw": "data.raw"},
        )
        method_b_node = ResolvedNode(
            id="method_b-MB-.default-D1",
            stage_id="method_b",
            module_id="MB",
            param_id=".default",
            module=FAKE_MODULE,
            inputs={"data_raw": "data/D1/.default/D1_data.csv"},
            outputs=["method_b/MB/.default/D1_result.csv"],
            input_name_mapping={"data_raw": "data.raw"},
        )
        gather_node = ResolvedNode(
            id="gather_aggregate-AGG-.default",
            stage_id="aggregate",
            module_id="AGG",
            param_id=".default",
            module=FAKE_MODULE,
            inputs={
                "input_0": "method_a/MA/.default/D1_result.csv",
                "input_1": "method_b/MB/.default/D1_result.csv",
            },
            outputs=["aggregate/AGG/.default/summary.csv"],
            input_name_mapping={
                "input_0": "method_output",
                "input_1": "method_output",
            },
        )

        all_nodes = [data_node, method_a_node, method_b_node, gather_node]

        gen = SnakemakeGenerator("GatherIntegrationBenchmark", "1.0", "Test")

        import tempfile

        with tempfile.NamedTemporaryFile(mode="w", suffix=".smk", delete=False) as f:
            tmp_path = Path(f.name)

        try:
            gen.generate_snakefile(
                nodes=all_nodes,
                collectors=[],
                output_path=tmp_path,
                debug_mode=True,
            )
            content = tmp_path.read_text()

            # Header
            assert "GatherIntegrationBenchmark" in content

            # All 4 rules
            assert "rule data_D1__default:" in content
            assert "rule method_a_MA__default_D1:" in content
            assert "rule method_b_MB__default_D1:" in content
            assert "rule gather_aggregate_AGG__default:" in content

            # 'all' rule references all outputs
            assert "rule all:" in content
            assert "data/D1/.default/D1_data.csv" in content
            assert "method_a/MA/.default/D1_result.csv" in content
            assert "method_b/MB/.default/D1_result.csv" in content
            assert "aggregate/AGG/.default/summary.csv" in content

            # Gather-specific: no benchmark directive, correct type
            # Find the gather rule section
            gather_idx = content.index("rule gather_aggregate")
            # Find the next rule after gather (the 'all' rule)
            all_idx = content.index("rule all:")
            gather_section = content[gather_idx:all_idx]

            assert "benchmark:" not in gather_section
            assert "# Type: Gather Stage" in content
        finally:
            tmp_path.unlink()


class TestGatherWithMultipleDatasets:
    """Integration test with multiple datasets flowing through providers."""

    def test_gather_collects_from_all_datasets_and_methods(self):
        """Gather should collect outputs from ALL (dataset x method) combinations."""
        bm_yaml = textwrap.dedent("""\
            id: MultiDatasetGather
            description: test
            version: "1.0"
            benchmarker: "Test"
            benchmark_yaml_spec: "0.3"
            software_backend: host
            software_environments:
              host:
                description: "Host"
            stages:
              - id: data
                modules:
                  - id: D1
                    software_environment: "host"
                    repository:
                      url: https://github.com/test/d
                      commit: aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
                  - id: D2
                    software_environment: "host"
                    repository:
                      url: https://github.com/test/d
                      commit: aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
                outputs:
                  - id: data.raw
                    path: "{dataset}_data.csv"
              - id: methods
                provides:
                  method_output: methods.result
                modules:
                  - id: MA
                    software_environment: "host"
                    repository:
                      url: https://github.com/test/m
                      commit: bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb
                  - id: MB
                    software_environment: "host"
                    repository:
                      url: https://github.com/test/m
                      commit: bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb
                inputs:
                  - data.raw
                outputs:
                  - id: methods.result
                    path: "{dataset}_result.csv"
              - id: gather_all
                modules:
                  - id: COLL
                    software_environment: "host"
                    repository:
                      url: https://github.com/test/c
                      commit: cccccccccccccccccccccccccccccccccccccccc
                inputs:
                  - gather: method_output
                outputs:
                  - id: gather_all.summary
                    path: summary.csv
        """)

        bm = Benchmark.from_yaml(bm_yaml)

        # Validate ordering
        _validate_gather_ordering(bm.stages)

        # Simulate 4 resolved nodes: 2 datasets x 2 methods
        method_nodes = []
        output_to_nodes = {}
        for dataset in ["D1", "D2"]:
            for method in ["MA", "MB"]:
                node_id = f"methods-{method}-.default-{dataset}"
                path = f"methods/{method}/.default/{dataset}_result.csv"
                method_nodes.append(
                    ResolvedNode(
                        id=node_id,
                        stage_id="methods",
                        module_id=method,
                        param_id=".default",
                        module=FAKE_MODULE,
                        outputs=[path],
                    )
                )
                output_to_nodes.setdefault("methods.result", []).append((node_id, path))

        inputs, name_mapping = _build_gather_inputs(
            ["method_output"], method_nodes, bm, output_to_nodes
        )

        # Should have 4 inputs (2 datasets x 2 methods)
        assert len(inputs) == 4
        assert all(v == "method_output" for v in name_mapping.values())

        # All 4 output paths should be present
        paths = set(inputs.values())
        assert "methods/MA/.default/D1_result.csv" in paths
        assert "methods/MA/.default/D2_result.csv" in paths
        assert "methods/MB/.default/D1_result.csv" in paths
        assert "methods/MB/.default/D2_result.csv" in paths


class TestGatherWithMultipleLabels:
    """Integration test for a gather stage that gathers from multiple labels."""

    def test_gather_two_labels(self):
        """A gather stage can collect from two different provides labels."""
        bm_yaml = textwrap.dedent("""\
            id: MultiLabelGather
            description: test
            version: "1.0"
            benchmarker: "Test"
            benchmark_yaml_spec: "0.3"
            software_backend: host
            software_environments:
              host:
                description: "Host"
            stages:
              - id: preprocessing
                provides:
                  preprocessed: preprocessing.output
                modules:
                  - id: PP
                    software_environment: "host"
                    repository:
                      url: https://github.com/test/pp
                      commit: aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
                outputs:
                  - id: preprocessing.output
                    path: preprocessed.csv
              - id: methods
                provides:
                  method_result: methods.result
                modules:
                  - id: MA
                    software_environment: "host"
                    repository:
                      url: https://github.com/test/ma
                      commit: bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb
                inputs:
                  - preprocessing.output
                outputs:
                  - id: methods.result
                    path: result.csv
              - id: summary
                modules:
                  - id: SUM
                    software_environment: "host"
                    repository:
                      url: https://github.com/test/sum
                      commit: cccccccccccccccccccccccccccccccccccccccc
                inputs:
                  - gather: preprocessed
                  - gather: method_result
                outputs:
                  - id: summary.out
                    path: summary.csv
        """)

        bm = Benchmark.from_yaml(bm_yaml)

        # Validate
        _validate_gather_ordering(bm.stages)

        # Verify multiple gather labels
        summary = bm.stages[2]
        assert summary.is_gather_stage()
        assert set(summary.get_gather_labels()) == {"preprocessed", "method_result"}

        # Simulate resolved nodes
        pp_node = ResolvedNode(
            id="preprocessing-PP-.default",
            stage_id="preprocessing",
            module_id="PP",
            param_id=".default",
            module=FAKE_MODULE,
            outputs=["preprocessing/PP/.default/preprocessed.csv"],
        )
        ma_node = ResolvedNode(
            id="methods-MA-.default",
            stage_id="methods",
            module_id="MA",
            param_id=".default",
            module=FAKE_MODULE,
            outputs=["methods/MA/.default/result.csv"],
        )

        output_to_nodes = {
            "preprocessing.output": [
                (
                    "preprocessing-PP-.default",
                    "preprocessing/PP/.default/preprocessed.csv",
                )
            ],
            "methods.result": [
                ("methods-MA-.default", "methods/MA/.default/result.csv")
            ],
        }

        inputs, name_mapping = _build_gather_inputs(
            ["preprocessed", "method_result"],
            [pp_node, ma_node],
            bm,
            output_to_nodes,
        )

        assert len(inputs) == 2
        assert inputs["input_0"] == "preprocessing/PP/.default/preprocessed.csv"
        assert inputs["input_1"] == "methods/MA/.default/result.csv"
        assert name_mapping["input_0"] == "preprocessed"
        assert name_mapping["input_1"] == "method_result"
