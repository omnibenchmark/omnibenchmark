"""
Unit tests for SnakemakeGenerator (omnibenchmark/backend/snakemake.py).

Covers:
- _sanitize_rule_name
- _make_human_name
- _write_header
- _write_all_rule
- _write_environment_directive (all backends)
- _write_resources_directive (no resources / with resources)
- _write_node_rule (regular, collector, gather; debug + exec modes)
- generate_snakefile (integration: writes a parseable file)
"""

import io
import json
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from omnibenchmark.backend._manifest import write_run_manifest
from omnibenchmark.backend.snakemake import SnakemakeGenerator, _make_human_name
from omnibenchmark.backend._snakemake_debug import DebugSnakemakeGenerator
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
    is_gather=False,
    is_collector=False,
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
        is_gather=is_gather,
        is_collector=is_collector,
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
    from omnibenchmark.model.params import Params

    return Params(d or {"method": "pca", "k": "10"})


def _gen(api_version=None) -> SnakemakeGenerator:
    kwargs = {
        "benchmark_name": "test-bench",
        "benchmark_version": "1.0",
        "benchmark_author": "tester",
    }
    if api_version is not None:
        kwargs["api_version"] = api_version
    return SnakemakeGenerator(**kwargs)


def _debug_gen() -> DebugSnakemakeGenerator:
    return DebugSnakemakeGenerator(
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
        name = _make_human_name(params)
        assert "k-10" in name or "method-pca" in name

    def test_unsafe_chars_replaced(self):
        params = _make_params({"path": "a/b", "q": "x y"})
        name = _make_human_name(params)
        assert "/" not in name
        assert " " not in name

    def test_long_name_truncated(self):
        # Build a param with a very long value
        params = _make_params({"k": "x" * 300})
        name = _make_human_name(params)
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
        out = _capture(_gen()._write_all_rule, [node])
        assert "rule all:" in out
        assert '"out/a.csv"' in out
        assert '"out/b.csv"' in out
        assert "default_target: True" in out

    def test_collector_node_outputs_included(self):
        node = _make_node(outputs=["out/result.csv"])
        collector_node = _make_node(outputs=["out/metrics.csv"], is_collector=True)
        out = _capture(_gen()._write_all_rule, [node, collector_node])
        assert '"out/metrics.csv"' in out

    def test_empty_nodes(self):
        out = _capture(_gen()._write_all_rule, [])
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
        out = _capture(_gen()._write_environment_directive, node.module)
        assert out == ""

    def test_conda(self):
        node = self._node_with_env("conda", "envs/myenv.yaml")
        out = _capture(_gen()._write_environment_directive, node.module)
        assert 'conda: "envs/myenv.yaml"' in out

    def test_apptainer(self):
        node = self._node_with_env("apptainer", "oras://my.registry/image:latest")
        out = _capture(_gen()._write_environment_directive, node.module)
        assert 'container: "oras://my.registry/image:latest"' in out

    def test_docker(self):
        node = self._node_with_env("docker", "docker://ubuntu:22.04")
        out = _capture(_gen()._write_environment_directive, node.module)
        assert 'container: "docker://ubuntu:22.04"' in out

    def test_envmodules(self):
        node = self._node_with_env("envmodules", "Python/3.10.0")
        out = _capture(_gen()._write_environment_directive, node.module)
        assert 'envmodules: "Python/3.10.0"' in out

    def test_host_backend_emits_nothing(self):
        node = self._node_with_env("host", "")
        out = _capture(_gen()._write_environment_directive, node.module)
        # host backend should produce no directive
        assert "conda:" not in out
        assert "container:" not in out
        assert "envmodules:" not in out


# ---------------------------------------------------------------------------
# _write_threads_directive
# ---------------------------------------------------------------------------


class TestWriteThreadsDirective:
    def test_no_resources_emits_nothing(self):
        """Default must stay unset: emitting threads would change scheduling
        for every existing benchmark."""
        node = _make_node(resources=None)
        out = _capture(_gen()._write_threads_directive, node)
        assert out == ""

    def test_cores_emits_threads(self):
        resources = MagicMock()
        resources.cores = 8
        node = _make_node(resources=resources)
        out = _capture(_gen()._write_threads_directive, node)
        assert "threads: 8" in out

    def test_zero_cores_emits_nothing(self):
        resources = MagicMock()
        resources.cores = 0
        node = _make_node(resources=resources)
        out = _capture(_gen()._write_threads_directive, node)
        assert out == ""


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
        out = _capture(_debug_gen()._write_node_rule, node)
        assert "rule stage1_mod1_default" in out
        assert "input:" in out
        assert "output:" in out
        assert "params:" in out
        assert "shell:" in out
        assert "touch" in out  # debug mode touches outputs

    def test_benchmark_directive_present_for_regular_node(self):
        node = _make_node(outputs=["out/stage1/mod1/default/result.csv"])
        out = _capture(_debug_gen()._write_node_rule, node)
        assert "benchmark:" in out

    def test_benchmark_directive_absent_for_collector(self):
        node = _make_node(is_collector=True)
        out = _capture(_debug_gen()._write_node_rule, node)
        assert "benchmark:" not in out

    def test_benchmark_directive_absent_for_gather(self):
        node = _make_node(is_gather=True)
        out = _capture(_debug_gen()._write_node_rule, node)
        assert "benchmark:" not in out

    def test_log_directive_present(self):
        node = _make_node()
        out = _capture(_debug_gen()._write_node_rule, node)
        assert ".logs/" in out
        assert ".log" in out

    def test_with_parameters(self):
        params = _make_params({"k": "10"})
        node = _make_node(parameters=params)
        out = _capture(_debug_gen()._write_node_rule, node)
        assert "cli_args=" in out
        # In debug mode a symlink line should appear
        assert "ln -sfn" in out

    def test_with_template_context(self):
        ctx = TemplateContext(provides={"dataset": "pbmc3k"})
        node = _make_node(template_context=ctx)
        out = _capture(_debug_gen()._write_node_rule, node)
        # --name is the module's own ID, not the dataset label
        assert "--name mod1" in out


# ---------------------------------------------------------------------------
# _write_node_rule — regular node, exec mode
# ---------------------------------------------------------------------------


class TestWriteNodeRuleExec:
    def test_no_touch_in_exec_mode(self):
        node = _make_node(outputs=["out/stage1/mod1/default/result.csv"])
        out = _capture(_gen()._write_node_rule, node)
        assert "touch" not in out

    def test_has_shebang_direct_exec(self):
        module = _make_module(has_shebang=True)
        node = _make_node(module=module, outputs=["out/stage1/mod1/default/result.csv"])
        out = _capture(_gen()._write_node_rule, node)
        assert "./{params.entrypoint}" in out

    def test_no_shebang_uses_interpreter(self):
        module = _make_module(has_shebang=False, interpreter="Rscript")
        node = _make_node(module=module, outputs=["out/stage1/mod1/default/result.csv"])
        out = _capture(_gen()._write_node_rule, node)
        assert "Rscript" in out

    def test_no_shebang_no_interpreter_defaults_python3(self):
        module = _make_module(has_shebang=False, interpreter=None)
        node = _make_node(module=module, outputs=["out/stage1/mod1/default/result.csv"])
        out = _capture(_gen()._write_node_rule, node)
        assert "python3" in out

    def test_inputs_resolved_to_absolute(self):
        node = _make_node(
            inputs={"data": "out/prev/data.csv"},
            outputs=["out/stage1/mod1/default/result.csv"],
            input_name_mapping={"data": "data"},
        )
        out = _capture(_gen()._write_node_rule, node)
        assert "INPUT_data=" in out
        assert "--data $INPUT_data" in out

    def test_parameters_json_written(self):
        params = _make_params({"k": "5"})
        node = _make_node(
            parameters=params,
            outputs=["out/stage1/mod1/default/result.csv"],
        )
        out = _capture(_gen()._write_node_rule, node)
        assert "parameters.json" in out

    def test_name_arg_is_module_id(self):
        # --name is the module ID for API >= 0.5 (spec §3.5)
        from omnibenchmark.model.benchmark import APIVersion

        node = _make_node(
            module_id="mymod",
            outputs=["out/stage1/mymod/default/result.csv"],
        )
        out = _capture(_gen(api_version=APIVersion.V0_5_0)._write_node_rule, node)
        assert "--name mymod" in out

    def test_name_arg_backward_compat_v0_4(self):
        # BUG WORKAROUND: For API <= 0.4, --name uses dataset name from path
        # This maintains backward compatibility with legacy modules that use
        # --name to construct output filenames ({name}_data.json)
        # Path structure: {prefix}/{stage}/{dataset}/{param}/...
        # parts[1] extracts the stage name in this case
        from omnibenchmark.model.benchmark import APIVersion

        node = _make_node(
            module_id="mymod",
            outputs=["out/stage1/dataset1/default/result.csv"],
        )
        out = _capture(_gen(api_version=APIVersion.V0_4_0)._write_node_rule, node)
        # With API 0.4, extracts parts[1] which is "stage1" in this path
        assert "--name stage1" in out
        assert "--name mymod" not in out

    def test_tee_logging(self):
        node = _make_node(outputs=["out/stage1/mod1/default/result.csv"])
        out = _capture(_gen()._write_node_rule, node)
        assert "tee" in out


# ---------------------------------------------------------------------------
# Collector/gather nodes written via _write_node_rule
# ---------------------------------------------------------------------------


class TestWriteNodeRuleCollectorV2:
    def _gather_node(self, **kwargs):
        return _make_node(
            node_id="gather_stage2-mod2-default",
            stage_id="stage2",
            inputs={"score_0": "out/a/score.csv", "score_1": "out/b/score.csv"},
            outputs=["out/gather/result.csv"],
            input_name_mapping={"score_0": "score", "score_1": "score"},
            is_gather=True,
            **kwargs,
        )

    def test_gather_debug_uses_v2(self):
        node = self._gather_node()
        out = _capture(_debug_gen()._write_node_rule, node)
        assert "GATHER STAGE" in out

    def test_gather_exec_uses_v2(self):
        node = self._gather_node()
        out = _capture(_gen()._write_node_rule, node)
        # v2 exec sets MODULE_DIR and OUTPUT_DIR variables
        assert "MODULE_DIR" in out
        assert "OUTPUT_DIR" in out

    def test_gather_exec_groups_inputs_by_name(self):
        node = self._gather_node()
        out = _capture(_gen()._write_node_rule, node)
        # Both inputs share original name 'score', so --score should appear once
        assert "--score" in out

    def test_gather_exec_shebang(self):
        module = _make_module(has_shebang=True)
        node = self._gather_node(module=module)
        out = _capture(_gen()._write_node_rule, node)
        assert "$MODULE_DIR/{params.entrypoint}" in out


# ---------------------------------------------------------------------------
# generate_snakefile (integration)
# ---------------------------------------------------------------------------


class TestGenerateSnakefile:
    def test_writes_file(self, tmp_path):
        gen = _debug_gen()
        node = _make_node(outputs=["out/stage1/mod1/default/result.csv"])
        output_file = tmp_path / "Snakefile"
        gen.generate_snakefile([node], output_file)
        assert output_file.exists()
        content = output_file.read_text()
        assert "rule all:" in content
        assert "rule stage1_mod1_default:" in content

    def test_collector_node_included(self, tmp_path):
        gen = _debug_gen()
        node = _make_node(outputs=["out/stage1/mod1/default/result.csv"])
        collector_node = _make_node(
            node_id="_collector_metrics-mod1-default",
            stage_id="_collector_metrics",
            outputs=["out/metrics/result.csv"],
            is_collector=True,
        )
        output_file = tmp_path / "Snakefile"
        gen.generate_snakefile([node, collector_node], output_file)
        content = output_file.read_text()
        assert "_collector_metrics" in content

    def test_exec_mode_no_touch(self, tmp_path):
        gen = _gen()
        node = _make_node(outputs=["out/stage1/mod1/default/result.csv"])
        output_file = tmp_path / "Snakefile"
        gen.generate_snakefile([node], output_file)
        content = output_file.read_text()
        assert "touch" not in content

    def test_empty_nodes(self, tmp_path):
        gen = _debug_gen()
        output_file = tmp_path / "Snakefile"
        gen.generate_snakefile([], output_file)
        content = output_file.read_text()
        assert "rule all:" in content


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


# ---------------------------------------------------------------------------
# Entrypoint invocation contract
#
# A shared repo checkout is made executable for every declared entrypoint by
# the resolver (chmod-all). Invocation then prefers the entrypoint's own
# shebang (direct ./exec) and falls back to the inferred interpreter only when
# there is no shebang.
#
# Two backend paths emit invocation commands: _write_shell (regular nodes) and
# _write_gather_shell (the shared aggregate path). Today only collectors
# (is_collector) reach the aggregate path in production; is_gather is dormant
# plumbing that routes to the same code and is not yet produced anywhere. The
# is_gather cases below pin the contract on that shared path so it can't drift
# from the regular path before gather is wired up by the explicit-gather
# feature coming after 0.7
# ---------------------------------------------------------------------------


@pytest.mark.short
class TestEntrypointInvocationContract:
    def test_regular_shebang_executes_directly(self):
        module = _make_module(entrypoint="select.R", has_shebang=True)
        node = _make_node(module=module, outputs=["out/s/m/default/r.csv"])
        out = _capture(_gen()._write_node_rule, node)
        assert "./{params.entrypoint}" in out
        assert "Rscript {params.entrypoint}" not in out

    def test_regular_no_shebang_falls_back_to_interpreter(self):
        module = _make_module(
            entrypoint="select.R", has_shebang=False, interpreter="Rscript"
        )
        node = _make_node(module=module, outputs=["out/s/m/default/r.csv"])
        out = _capture(_gen()._write_node_rule, node)
        assert "Rscript {params.entrypoint}" in out
        assert "./{params.entrypoint}" not in out

    def test_regular_no_shebang_no_interpreter_defaults_python3(self):
        module = _make_module(entrypoint="run", has_shebang=False, interpreter=None)
        node = _make_node(module=module, outputs=["out/s/m/default/r.csv"])
        out = _capture(_gen()._write_node_rule, node)
        assert "python3 {params.entrypoint}" in out

    def test_gather_shebang_executes_directly(self):
        module = _make_module(entrypoint="integrate.R", has_shebang=True)
        node = _make_node(
            node_id="gather_s2-m2-default",
            stage_id="s2",
            module=module,
            inputs={"score_0": "out/a/score.csv"},
            outputs=["out/gather/result.csv"],
            input_name_mapping={"score_0": "score"},
            is_gather=True,
        )
        out = _capture(_gen()._write_node_rule, node)
        assert "$MODULE_DIR/{params.entrypoint}" in out
        assert "Rscript $MODULE_DIR" not in out

    def test_gather_no_shebang_falls_back_to_interpreter(self):
        module = _make_module(
            entrypoint="integrate.R", has_shebang=False, interpreter="Rscript"
        )
        node = _make_node(
            node_id="gather_s2-m2-default",
            stage_id="s2",
            module=module,
            inputs={"score_0": "out/a/score.csv"},
            outputs=["out/gather/result.csv"],
            input_name_mapping={"score_0": "score"},
            is_gather=True,
        )
        out = _capture(_gen()._write_node_rule, node)
        assert "Rscript $MODULE_DIR/{params.entrypoint}" in out


# ---------------------------------------------------------------------------
# lineage.json — provenance for nodes whose path cannot encode their ancestry
# (design 010 §3.3/§5.2).
# ---------------------------------------------------------------------------


def _gen_with(nodes):
    gen = SnakemakeGenerator("b", "1.0", "me")
    gen._nodes_by_id = {n.id: n for n in nodes}
    return gen


@pytest.mark.short
def test_no_lineage_sidecar_for_linear_nodes():
    """A linear node's directory prefix already IS its lineage — writing a
    sidecar there would restate the path for every node in the plan."""
    node = _make_node(node_id="s-m-default")
    gen = _gen_with([node])

    assert gen._lineage_record(node) is None
    assert gen._lineage_sidecar_lines(node) == []


@pytest.mark.short
def test_join_lineage_sidecar_records_every_branch():
    """A join's path follows only its deepest input, so the other branches are
    recoverable from the filesystem only if written down."""
    deep = _make_node(
        node_id="root-deep",
        stage_id="deep",
        module_id="D1",
        outputs=["root/R1/deep/D1/.h1/d.json"],
    )
    shallow = _make_node(
        node_id="root-shallow",
        stage_id="shallow",
        module_id="S1",
        outputs=["root/R1/shallow/S1/.h2/s.json"],
    )
    join = _make_node(node_id="join-J1-abcd1234", stage_id="join", module_id="J1")
    object.__setattr__(join, "parents", [deep.id, shallow.id])

    record = _gen_with([deep, shallow, join])._lineage_record(join)

    assert record["kind"] == "join"
    assert [m["stage"] for m in record["members"]] == ["deep", "shallow"]
    # the branch the path drops must be present, with enough to locate it
    shallow_member = record["members"][1]
    assert shallow_member["module"] == "S1"
    assert shallow_member["dir"] == "root/R1/shallow/S1/.h2"
    assert shallow_member["commit"] == "abc1234"


@pytest.mark.short
def test_gather_lineage_sidecar_uses_gathered_from():
    """A gather cuts the chain entirely: its members appear nowhere in the
    path, and `gathered_from` is the only record of them."""
    m1 = _make_node(node_id="m1", stage_id="methods", module_id="M1")
    m2 = _make_node(node_id="m2", stage_id="methods", module_id="M2")
    g = _make_node(node_id="agg-G1", stage_id="agg", module_id="G1", is_gather=True)
    object.__setattr__(g, "gathered_from", [m1.id, m2.id])

    record = _gen_with([m1, m2, g])._lineage_record(g)

    assert record["kind"] == "gather"
    assert [m["module"] for m in record["members"]] == ["M1", "M2"]


@pytest.mark.short
def test_lineage_sidecar_is_one_shell_line_with_escaped_braces():
    """Every line of a shell block is indented, so an indented heredoc
    terminator would not close the heredoc — the record goes out as a single
    echo, with braces doubled for Snakemake's own formatting pass."""
    a = _make_node(node_id="a", stage_id="a", module_id="A")
    b = _make_node(node_id="b", stage_id="b", module_id="B")
    join = _make_node(node_id="j", stage_id="j", module_id="J")
    object.__setattr__(join, "parents", [a.id, b.id])

    lines = _gen_with([a, b, join])._lineage_sidecar_lines(join)

    assert len(lines) == 1
    assert lines[0].endswith("> $OUTPUT_DIR/lineage.json")
    assert "{{" in lines[0] and "}}" in lines[0]
    # round-trips once Snakemake collapses the doubled braces
    payload = lines[0].split("echo '", 1)[1].rsplit("' >", 1)[0]
    assert json.loads(payload.replace("{{", "{").replace("}}", "}"))["kind"] == "join"


@pytest.mark.short
def test_metadata_records_every_module_not_just_every_repo(tmp_path):
    """One repository commonly hosts several modules, distinguished only by
    entrypoint (scrapper -> filter/normalize/select/pca). Deduping the record
    on (repo, commit) dropped all but the first, so modules.txt under-reported
    what actually ran."""
    from omnibenchmark.backend._metadata import save_metadata

    shared = dict(repo_url="https://github.com/x/scrapper", commit="18843a7")
    nodes = [
        _make_node(
            node_id="filt-fi",
            stage_id="FILT",
            module_id="fi-scrapper",
            module=_make_module(entrypoint="filter", **shared),
        ),
        _make_node(
            node_id="norm-nr",
            stage_id="NORM",
            module_id="nr-scrapper",
            module=_make_module(entrypoint="normalize", **shared),
        ),
    ]
    yaml_path = tmp_path / "benchmark.yaml"
    yaml_path.write_text("id: b\n")

    save_metadata(benchmark_yaml_path=yaml_path, output_dir=tmp_path, nodes=nodes)
    listed = (tmp_path / ".metadata" / "modules.txt").read_text()

    assert "FILT/fi-scrapper" in listed
    assert "NORM/nr-scrapper" in listed
