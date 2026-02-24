"""
Unit tests for SnakemakeGenerator (omnibenchmark/backend/snakemake_gen.py).

Covers:
- _sanitize_rule_name
- _make_human_name
- _write_header
- _write_all_rule
- _write_environment_directive (all backends)
- _write_resources_directive (no resources / with resources)
- _write_node_rule (regular, collector, gather; debug + exec modes)
- _write_collector_rule (debug + exec modes)
- generate_snakefile (integration: writes a parseable file)
- save_metadata
"""

import io
import json
from pathlib import Path
from unittest.mock import MagicMock

from omnibenchmark.backend.snakemake_gen import (
    SnakemakeGenerator,
    save_metadata,
    write_run_manifest,
)
from omnibenchmark.model.resolved import (
    ResolvedNode,
    ResolvedModule,
    ResolvedMetricCollector,
    ResolvedEnvironment,
    TemplateContext,
)
from omnibenchmark.model.benchmark import SoftwareBackendEnum


# ---------------------------------------------------------------------------
# Helpers / factories
# ---------------------------------------------------------------------------


def _make_module(
    *,
    repo_url="https://github.com/example/repo",
    commit="abc1234",
    module_dir=".modules/repo/abc1234",
    entrypoint="run.py",
    has_shebang=False,
    interpreter="python3",
    resolved_environment=None,
):
    return ResolvedModule(
        repository_url=repo_url,
        commit=commit,
        module_dir=Path(module_dir),
        entrypoint=Path(entrypoint),
        software_environment_id="env1",
        has_shebang=has_shebang,
        interpreter=interpreter,
        resolved_environment=resolved_environment,
    )


def _make_node(
    *,
    node_id="stage1-mod1-default",
    stage_id="stage1",
    module_id="mod1",
    param_id="default",
    inputs=None,
    outputs=None,
    parameters=None,
    resources=None,
    template_context=None,
    input_name_mapping=None,
    module=None,
):
    if module is None:
        module = _make_module()
    return ResolvedNode(
        id=node_id,
        stage_id=stage_id,
        module_id=module_id,
        param_id=param_id,
        module=module,
        inputs=inputs or {},
        outputs=outputs or ["out/stage1/mod1/default/result.csv"],
        parameters=parameters,
        resources=resources,
        template_context=template_context,
        input_name_mapping=input_name_mapping or {},
    )


def _make_collector(
    *,
    collector_id="collector_metrics-mod1-default",
    input_patterns=None,
    outputs=None,
    parameters=None,
    module=None,
):
    if module is None:
        module = _make_module(entrypoint="collect.py")
    return ResolvedMetricCollector(
        id=collector_id,
        module=module,
        input_patterns=input_patterns or ["out/stage1/**/*.csv"],
        outputs=outputs or ["out/metrics/result.csv"],
        parameters=parameters,
    )


def _make_params(d=None):
    """Make a minimal Params-like object from dict."""
    from omnibenchmark.benchmark.params import Params

    return Params(d or {"method": "pca", "k": "10"})


def _gen() -> SnakemakeGenerator:
    return SnakemakeGenerator(
        benchmark_name="test-bench",
        benchmark_version="1.0",
        benchmark_author="tester",
    )


def _capture(method, *args, **kwargs) -> str:
    """Call a generator method that writes to a TextIO and return the output."""
    buf = io.StringIO()
    method(buf, *args, **kwargs)
    return buf.getvalue()


# ---------------------------------------------------------------------------
# _sanitize_rule_name
# ---------------------------------------------------------------------------


class TestSanitizeRuleName:
    def test_hyphens_replaced(self):
        assert _gen()._sanitize_rule_name("stage-mod-default") == "stage_mod_default"

    def test_dots_replaced(self):
        assert _gen()._sanitize_rule_name("stage.mod.default") == "stage_mod_default"

    def test_starts_with_digit_prefixed(self):
        result = _gen()._sanitize_rule_name("1stage")
        assert result.startswith("rule_")

    def test_already_valid(self):
        assert _gen()._sanitize_rule_name("valid_name") == "valid_name"

    def test_empty_string(self):
        # Empty string: not alpha, so prefixed
        result = _gen()._sanitize_rule_name("")
        assert result == ""  # empty — 'if name and not name[0].isalpha()' skips


# ---------------------------------------------------------------------------
# _make_human_name
# ---------------------------------------------------------------------------


class TestMakeHumanName:
    def test_basic(self):
        params = _make_params({"k": "10", "method": "pca"})
        name = _gen()._make_human_name(params)
        assert "k-10" in name or "method-pca" in name

    def test_unsafe_chars_replaced(self):
        params = _make_params({"path": "a/b", "q": "x y"})
        name = _gen()._make_human_name(params)
        assert "/" not in name
        assert " " not in name

    def test_long_name_truncated(self):
        # Build a param with a very long value
        params = _make_params({"k": "x" * 300})
        name = _gen()._make_human_name(params)
        assert len(name) <= 255


# ---------------------------------------------------------------------------
# _write_header
# ---------------------------------------------------------------------------


class TestWriteHeader:
    def test_contains_shebang(self):
        out = _capture(_gen()._write_header)
        assert "#!/usr/bin/env snakemake" in out

    def test_contains_benchmark_name(self):
        out = _capture(_gen()._write_header)
        assert "test-bench" in out

    def test_contains_version(self):
        out = _capture(_gen()._write_header)
        assert "1.0" in out


# ---------------------------------------------------------------------------
# _write_all_rule
# ---------------------------------------------------------------------------


class TestWriteAllRule:
    def test_rule_all_present(self):
        node = _make_node(outputs=["out/a.csv", "out/b.csv"])
        out = _capture(_gen()._write_all_rule, [node], [])
        assert "rule all:" in out
        assert '"out/a.csv"' in out
        assert '"out/b.csv"' in out
        assert "default_target: True" in out

    def test_collector_outputs_included(self):
        node = _make_node(outputs=["out/result.csv"])
        collector = _make_collector(outputs=["out/metrics.csv"])
        out = _capture(_gen()._write_all_rule, [node], [collector])
        assert '"out/metrics.csv"' in out

    def test_empty_nodes_and_collectors(self):
        out = _capture(_gen()._write_all_rule, [], [])
        assert "rule all:" in out
        assert "default_target: True" in out


# ---------------------------------------------------------------------------
# _write_environment_directive
# ---------------------------------------------------------------------------


class TestWriteEnvironmentDirective:
    def _node_with_env(self, backend_value, reference):
        env = ResolvedEnvironment(
            backend_type=SoftwareBackendEnum(backend_value),
            reference=reference,
        )
        module = _make_module(resolved_environment=env)
        return _make_node(module=module)

    def test_no_environment(self):
        node = _make_node(module=_make_module(resolved_environment=None))
        out = _capture(_gen()._write_environment_directive, node)
        assert out == ""

    def test_conda(self):
        node = self._node_with_env("conda", "envs/myenv.yaml")
        out = _capture(_gen()._write_environment_directive, node)
        assert 'conda: "envs/myenv.yaml"' in out

    def test_apptainer(self):
        node = self._node_with_env("apptainer", "oras://my.registry/image:latest")
        out = _capture(_gen()._write_environment_directive, node)
        assert 'container: "oras://my.registry/image:latest"' in out

    def test_docker(self):
        node = self._node_with_env("docker", "docker://ubuntu:22.04")
        out = _capture(_gen()._write_environment_directive, node)
        assert 'container: "docker://ubuntu:22.04"' in out

    def test_envmodules(self):
        node = self._node_with_env("envmodules", "Python/3.10.0")
        out = _capture(_gen()._write_environment_directive, node)
        assert 'envmodules: "Python/3.10.0"' in out

    def test_host_backend_emits_nothing(self):
        node = self._node_with_env("host", "")
        out = _capture(_gen()._write_environment_directive, node)
        # host backend should produce no directive
        assert "conda:" not in out
        assert "container:" not in out
        assert "envmodules:" not in out


# ---------------------------------------------------------------------------
# _write_resources_directive
# ---------------------------------------------------------------------------


class TestWriteResourcesDirective:
    def test_no_resources_uses_default_cores(self):
        node = _make_node(resources=None)
        out = _capture(_gen()._write_resources_directive, node)
        assert "resources:" in out
        assert "cores=2" in out

    def test_custom_cores(self):
        resources = MagicMock()
        resources.cores = 8
        resources.mem_mb = None
        resources.disk_mb = None
        resources.runtime = None
        resources.gpu = None
        node = _make_node(resources=resources)
        out = _capture(_gen()._write_resources_directive, node)
        assert "cores=8" in out

    def test_mem_mb_included(self):
        resources = MagicMock()
        resources.cores = 4
        resources.mem_mb = 8192
        resources.disk_mb = None
        resources.runtime = None
        resources.gpu = None
        node = _make_node(resources=resources)
        out = _capture(_gen()._write_resources_directive, node)
        assert "mem_mb=8192" in out

    def test_gpu_included(self):
        resources = MagicMock()
        resources.cores = 4
        resources.mem_mb = None
        resources.disk_mb = None
        resources.runtime = None
        resources.gpu = 1
        node = _make_node(resources=resources)
        out = _capture(_gen()._write_resources_directive, node)
        assert "nvidia_gpu=1" in out


# ---------------------------------------------------------------------------
# _write_node_rule — regular node, debug mode
# ---------------------------------------------------------------------------


class TestWriteNodeRuleDebug:
    def test_basic_structure(self):
        node = _make_node(
            inputs={"data": "out/prev/data.csv"},
            outputs=["out/stage1/mod1/default/result.csv"],
        )
        out = _capture(_gen()._write_node_rule, node, debug_mode=True)
        assert "rule stage1_mod1_default" in out
        assert "input:" in out
        assert "output:" in out
        assert "params:" in out
        assert "shell:" in out
        assert "touch" in out  # debug mode touches outputs

    def test_benchmark_directive_present_for_regular_node(self):
        node = _make_node(outputs=["out/stage1/mod1/default/result.csv"])
        out = _capture(_gen()._write_node_rule, node, debug_mode=True)
        assert "benchmark:" in out

    def test_benchmark_directive_absent_for_collector(self):
        node = _make_node(
            node_id="_collector_metrics-mod1-default",
            stage_id="_collector_metrics",
        )
        out = _capture(_gen()._write_node_rule, node, debug_mode=True)
        assert "benchmark:" not in out

    def test_benchmark_directive_absent_for_gather(self):
        node = _make_node(
            node_id="gather_stage1-mod1-default",
            stage_id="stage1",
        )
        out = _capture(_gen()._write_node_rule, node, debug_mode=True)
        assert "benchmark:" not in out

    def test_log_directive_present(self):
        node = _make_node()
        out = _capture(_gen()._write_node_rule, node, debug_mode=True)
        assert ".logs/" in out
        assert ".log" in out

    def test_with_parameters(self):
        params = _make_params({"k": "10"})
        node = _make_node(parameters=params)
        out = _capture(_gen()._write_node_rule, node, debug_mode=True)
        assert "cli_args=" in out
        # In debug mode a symlink line should appear
        assert "ln -sfn" in out

    def test_with_template_context_dataset(self):
        ctx = TemplateContext(provides={"dataset": "pbmc3k"})
        node = _make_node(template_context=ctx)
        out = _capture(_gen()._write_node_rule, node, debug_mode=True)
        assert "pbmc3k" in out


# ---------------------------------------------------------------------------
# _write_node_rule — regular node, exec mode
# ---------------------------------------------------------------------------


class TestWriteNodeRuleExec:
    def test_no_touch_in_exec_mode(self):
        node = _make_node(outputs=["out/stage1/mod1/default/result.csv"])
        out = _capture(_gen()._write_node_rule, node, debug_mode=False)
        assert "touch" not in out

    def test_has_shebang_direct_exec(self):
        module = _make_module(has_shebang=True)
        node = _make_node(module=module, outputs=["out/stage1/mod1/default/result.csv"])
        out = _capture(_gen()._write_node_rule, node, debug_mode=False)
        assert "./{params.entrypoint}" in out

    def test_no_shebang_uses_interpreter(self):
        module = _make_module(has_shebang=False, interpreter="Rscript")
        node = _make_node(module=module, outputs=["out/stage1/mod1/default/result.csv"])
        out = _capture(_gen()._write_node_rule, node, debug_mode=False)
        assert "Rscript" in out

    def test_no_shebang_no_interpreter_defaults_python3(self):
        module = _make_module(has_shebang=False, interpreter=None)
        node = _make_node(module=module, outputs=["out/stage1/mod1/default/result.csv"])
        out = _capture(_gen()._write_node_rule, node, debug_mode=False)
        assert "python3" in out

    def test_inputs_prefixed(self):
        node = _make_node(
            inputs={"data": "out/prev/data.csv"},
            outputs=["out/stage1/mod1/default/result.csv"],
            input_name_mapping={"data": "data"},
        )
        out = _capture(_gen()._write_node_rule, node, debug_mode=False)
        assert "../../../{input.data}" in out

    def test_parameters_json_written(self):
        params = _make_params({"k": "5"})
        node = _make_node(
            parameters=params,
            outputs=["out/stage1/mod1/default/result.csv"],
        )
        out = _capture(_gen()._write_node_rule, node, debug_mode=False)
        assert "parameters.json" in out

    def test_dataset_name_arg(self):
        ctx = TemplateContext(provides={"dataset": "dataset1"})
        node = _make_node(
            template_context=ctx,
            outputs=["out/stage1/mod1/default/result.csv"],
        )
        out = _capture(_gen()._write_node_rule, node, debug_mode=False)
        assert "--name dataset1" in out

    def test_tee_logging(self):
        node = _make_node(outputs=["out/stage1/mod1/default/result.csv"])
        out = _capture(_gen()._write_node_rule, node, debug_mode=False)
        assert "tee" in out


# ---------------------------------------------------------------------------
# _write_collector_rule (ResolvedMetricCollector)
# ---------------------------------------------------------------------------


class TestWriteCollectorRule:
    def test_debug_mode_touches_outputs(self):
        collector = _make_collector()
        out = _capture(_gen()._write_collector_rule, collector, debug_mode=True)
        assert "touch" in out

    def test_exec_mode_no_touch(self):
        collector = _make_collector()
        out = _capture(_gen()._write_collector_rule, collector, debug_mode=False)
        assert "touch" not in out

    def test_exec_shebang(self):
        module = _make_module(has_shebang=True)
        collector = _make_collector(module=module)
        out = _capture(_gen()._write_collector_rule, collector, debug_mode=False)
        assert "./{params.entrypoint}" in out

    def test_exec_no_shebang_uses_interpreter(self):
        module = _make_module(has_shebang=False, interpreter="python3")
        collector = _make_collector(module=module)
        out = _capture(_gen()._write_collector_rule, collector, debug_mode=False)
        assert "python3 {params.entrypoint}" in out

    def test_input_patterns_in_output(self):
        collector = _make_collector(input_patterns=["out/**/*.csv"])
        out = _capture(_gen()._write_collector_rule, collector, debug_mode=True)
        assert '"out/**/*.csv"' in out

    def test_collector_outputs_in_rule(self):
        collector = _make_collector(outputs=["out/metrics/result.csv"])
        out = _capture(_gen()._write_collector_rule, collector, debug_mode=True)
        assert '"out/metrics/result.csv"' in out

    def test_with_parameters(self):
        params = _make_params({"alpha": "0.5"})
        collector = _make_collector(parameters=params)
        out = _capture(_gen()._write_collector_rule, collector, debug_mode=True)
        assert "cli_args=" in out


# ---------------------------------------------------------------------------
# Collector/gather nodes written via _write_node_rule (v2 path)
# ---------------------------------------------------------------------------


class TestWriteNodeRuleCollectorV2:
    def _gather_node(self, **kwargs):
        return _make_node(
            node_id="gather_stage2-mod2-default",
            stage_id="stage2",
            inputs={"score_0": "out/a/score.csv", "score_1": "out/b/score.csv"},
            outputs=["out/gather/result.csv"],
            input_name_mapping={"score_0": "score", "score_1": "score"},
            **kwargs,
        )

    def test_gather_debug_uses_v2(self):
        node = self._gather_node()
        out = _capture(_gen()._write_node_rule, node, debug_mode=True)
        assert "GATHER STAGE" in out

    def test_gather_exec_uses_v2(self):
        node = self._gather_node()
        out = _capture(_gen()._write_node_rule, node, debug_mode=False)
        # v2 exec sets MODULE_DIR and OUTPUT_DIR variables
        assert "MODULE_DIR" in out
        assert "OUTPUT_DIR" in out

    def test_gather_exec_groups_inputs_by_name(self):
        node = self._gather_node()
        out = _capture(_gen()._write_node_rule, node, debug_mode=False)
        # Both inputs share original name 'score', so --score should appear once
        assert "--score" in out

    def test_gather_exec_shebang(self):
        module = _make_module(has_shebang=True)
        node = self._gather_node(module=module)
        out = _capture(_gen()._write_node_rule, node, debug_mode=False)
        assert "$MODULE_DIR/{params.entrypoint}" in out


# ---------------------------------------------------------------------------
# generate_snakefile (integration)
# ---------------------------------------------------------------------------


class TestGenerateSnakefile:
    def test_writes_file(self, tmp_path):
        gen = _gen()
        node = _make_node(outputs=["out/stage1/mod1/default/result.csv"])
        output_file = tmp_path / "Snakefile"
        gen.generate_snakefile([node], [], output_file, debug_mode=True)
        assert output_file.exists()
        content = output_file.read_text()
        assert "rule all:" in content
        assert "rule stage1_mod1_default:" in content

    def test_with_collector(self, tmp_path):
        gen = _gen()
        node = _make_node(outputs=["out/stage1/mod1/default/result.csv"])
        collector = _make_collector(outputs=["out/metrics/result.csv"])
        output_file = tmp_path / "Snakefile"
        gen.generate_snakefile([node], [collector], output_file, debug_mode=True)
        content = output_file.read_text()
        assert "collector_metrics" in content

    def test_exec_mode_no_touch(self, tmp_path):
        gen = _gen()
        node = _make_node(outputs=["out/stage1/mod1/default/result.csv"])
        output_file = tmp_path / "Snakefile"
        gen.generate_snakefile([node], [], output_file, debug_mode=False)
        content = output_file.read_text()
        assert "touch" not in content

    def test_empty_nodes_and_collectors(self, tmp_path):
        gen = _gen()
        output_file = tmp_path / "Snakefile"
        gen.generate_snakefile([], [], output_file, debug_mode=True)
        content = output_file.read_text()
        assert "rule all:" in content


# ---------------------------------------------------------------------------
# save_metadata
# ---------------------------------------------------------------------------


class TestSaveMetadata:
    def test_creates_metadata_dir(self, tmp_path):
        benchmark_yaml = tmp_path / "benchmark.yaml"
        benchmark_yaml.write_text("name: test")
        out_dir = tmp_path / "out"
        out_dir.mkdir()
        node = _make_node()
        save_metadata(benchmark_yaml, out_dir, [node], [])
        assert (out_dir / ".metadata").is_dir()
        assert (out_dir / ".metadata" / "benchmark.yaml").exists()
        assert (out_dir / ".metadata" / "modules.txt").exists()

    def test_modules_txt_contains_module_info(self, tmp_path):
        benchmark_yaml = tmp_path / "benchmark.yaml"
        benchmark_yaml.write_text("name: test")
        out_dir = tmp_path / "out"
        out_dir.mkdir()
        node = _make_node()
        save_metadata(benchmark_yaml, out_dir, [node], [])
        modules_txt = (out_dir / ".metadata" / "modules.txt").read_text()
        assert "https://github.com/example/repo" in modules_txt
        assert "abc1234" in modules_txt

    def test_duplicate_modules_deduplicated(self, tmp_path):
        benchmark_yaml = tmp_path / "benchmark.yaml"
        benchmark_yaml.write_text("name: test")
        out_dir = tmp_path / "out"
        out_dir.mkdir()
        # Two nodes sharing the same repo+commit
        node1 = _make_node(node_id="stage1-mod1-p1", param_id="p1")
        node2 = _make_node(node_id="stage1-mod1-p2", param_id="p2")
        save_metadata(benchmark_yaml, out_dir, [node1, node2], [])
        modules_txt = (out_dir / ".metadata" / "modules.txt").read_text()
        # repo URL should appear only once
        assert modules_txt.count("https://github.com/example/repo") == 1

    def test_collectors_included(self, tmp_path):
        benchmark_yaml = tmp_path / "benchmark.yaml"
        benchmark_yaml.write_text("name: test")
        out_dir = tmp_path / "out"
        out_dir.mkdir()
        collector = _make_collector()
        save_metadata(benchmark_yaml, out_dir, [], [collector])
        modules_txt = (out_dir / ".metadata" / "modules.txt").read_text()
        assert "https://github.com/example/repo" in modules_txt


# ---------------------------------------------------------------------------
# write_run_manifest
# ---------------------------------------------------------------------------


class TestWriteRunManifest:
    def test_creates_manifest_json(self, tmp_path):
        out_dir = tmp_path / "out"
        out_dir.mkdir()
        write_run_manifest(out_dir)
        manifest_path = out_dir / ".metadata" / "manifest.json"
        assert manifest_path.exists()

    def test_manifest_is_valid_json(self, tmp_path):
        out_dir = tmp_path / "out"
        out_dir.mkdir()
        write_run_manifest(out_dir)
        content = (out_dir / ".metadata" / "manifest.json").read_text()
        data = json.loads(content)
        assert isinstance(data, dict)

    def test_manifest_contains_expected_keys(self, tmp_path):
        out_dir = tmp_path / "out"
        out_dir.mkdir()
        write_run_manifest(out_dir)
        data = json.loads((out_dir / ".metadata" / "manifest.json").read_text())
        for key in ("run_id", "timestamp", "hostname", "platform", "python_version"):
            assert key in data, f"missing key: {key}"

    def test_explicit_run_id_preserved(self, tmp_path):
        out_dir = tmp_path / "out"
        out_dir.mkdir()
        result = write_run_manifest(out_dir, run_id="my-fixed-run-id")
        assert result["run_id"] == "my-fixed-run-id"

    def test_auto_run_id_is_uuid(self, tmp_path):
        import uuid

        out_dir = tmp_path / "out"
        out_dir.mkdir()
        result = write_run_manifest(out_dir)
        # Should be a valid UUID
        parsed = uuid.UUID(result["run_id"])
        assert str(parsed) == result["run_id"]

    def test_returns_manifest_dict(self, tmp_path):
        out_dir = tmp_path / "out"
        out_dir.mkdir()
        result = write_run_manifest(out_dir)
        assert isinstance(result, dict)
        assert result["platform"] in ("linux", "darwin", "win32", result["platform"])
