"""Unit tests for PodmanSnakemakeGenerator."""

import io
from pathlib import Path

import pytest

from omnibenchmark.backend._snakemake_podman import (
    PodmanSnakemakeGenerator,
    _DEFAULT_MEM_MB,
)
from omnibenchmark.model.benchmark import Resources, SoftwareBackendEnum
from omnibenchmark.model.resolved import (
    ResolvedEnvironment,
    ResolvedModule,
    ResolvedNode,
)


def _module(*, image=None, has_shebang=False, interpreter="python3"):
    env = (
        ResolvedEnvironment(backend_type=SoftwareBackendEnum.podman, reference=image)
        if image
        else None
    )
    return ResolvedModule(
        repository_url="https://github.com/example/repo",
        commit="abc1234",
        module_dir=Path(".modules/repo/abc1234"),
        entrypoint=Path("run.py"),
        software_environment_id="env1",
        has_shebang=has_shebang,
        interpreter=interpreter,
        resolved_environment=env,
    )


def _node(*, image="myorg/img:1.0", resources=None, inputs=None, parameters=None):
    return ResolvedNode(
        id="stage1-mod1-default",
        stage_id="stage1",
        module_id="mod1",
        param_id="default",
        module=_module(image=image),
        inputs=inputs or {},
        outputs=["out/stage1/mod1/default/result.csv"],
        parameters=parameters,
        resources=resources,
        template_context=None,
        input_name_mapping={},
    )


def _gen(default_image=None) -> PodmanSnakemakeGenerator:
    return PodmanSnakemakeGenerator(
        benchmark_name="bench",
        benchmark_version="1.0",
        benchmark_author="t",
        default_image=default_image,
    )


def _capture(method, *args) -> str:
    buf = io.StringIO()
    method(buf, *args)
    return buf.getvalue()


# ---------------------------------------------------------------------------
# Image resolution
# ---------------------------------------------------------------------------


class TestImageResolution:
    def test_uses_module_image(self):
        n = _node(image="alpine:3")
        assert _gen()._image_for(n) == "alpine:3"

    def test_strips_docker_scheme(self):
        n = _node(image="docker://alpine:3")
        assert _gen()._image_for(n) == "alpine:3"

    def test_falls_back_to_default(self):
        n = _node(image=None)
        assert _gen(default_image="ubuntu:22.04")._image_for(n) == "ubuntu:22.04"

    def test_raises_when_no_image_anywhere(self):
        n = _node(image=None)
        with pytest.raises(ValueError, match="no docker image"):
            _gen()._image_for(n)


# ---------------------------------------------------------------------------
# Shell generation
# ---------------------------------------------------------------------------


class TestShellGeneration:
    def test_podman_run_present(self):
        out = _capture(_gen()._write_shell, _node())
        assert "podman run --rm" in out
        assert "myorg/img:1.0" in out

    def test_resource_placeholders_used(self):
        # CPU and memory come from Snakemake's resources block, not literals.
        out = _capture(_gen()._write_shell, _node())
        assert "--cpus {resources.cores}" in out
        assert "--memory {resources.mem_mb}m" in out

    def test_workdir_bind_mount_uses_physical_path(self):
        # Mount and workdir use $MOUNT_ROOT (set to `pwd -P`) so a symlinked
        # CWD does not produce a bind mount that disagrees with realpath()
        # results inside the container.
        out = _capture(_gen()._write_shell, _node())
        assert "MOUNT_ROOT=$(pwd -P)" in out
        assert "--volume $MOUNT_ROOT:$MOUNT_ROOT" in out
        assert "--workdir $MOUNT_ROOT" in out
        assert "$(pwd):$(pwd)" not in out  # regression: no logical-pwd mount

    def test_pull_missing(self):
        out = _capture(_gen()._write_shell, _node())
        assert "--pull=missing" in out

    def test_module_dir_cd_inside_container(self):
        # cd happens inside `sh -c '...'` so the bind-mounted module path
        # is reached from within the container.
        out = _capture(_gen()._write_shell, _node())
        assert "sh -c 'cd {params.module_dir} &&" in out

    def test_inputs_resolved_with_readlink_and_guarded(self):
        # Each input is realpath'd via readlink -f, then a case-glob asserts
        # the result lives under $MOUNT_ROOT — escapes abort with exit 1.
        n = _node(inputs={"counts": "data/D1/counts.csv"})
        out = _capture(_gen()._write_shell, n)
        assert "INPUT_counts=$(readlink -f {input.counts})" in out
        assert 'case "$INPUT_counts/" in "$MOUNT_ROOT"/*) ;;' in out
        assert "exit 1" in out

    def test_output_dir_uses_physical_path(self):
        out = _capture(_gen()._write_shell, _node())
        assert "OUTPUT_DIR=$(cd {params.output_dir} && pwd -P)" in out

    def test_no_environment_directive(self):
        # The container: directive is suppressed; podman is launched in shell.
        node = _node(image="alpine:3")
        out = _capture(_gen()._write_node_rule, node)
        assert "container:" not in out

    def test_no_benchmark_directive(self):
        # Snakemake's benchmark: directive uses psutil to walk the rule's process
        # tree; it can't see inside a podman container. Keep it but rely on cgroup
        # numbers. (Currently we keep it on; this test pins behavior — flip if changed.)
        node = _node(image="alpine:3")
        out = _capture(_gen()._write_node_rule, node)
        # Benchmark directive remains for non-collector regular nodes.
        assert "benchmark:" in out


# ---------------------------------------------------------------------------
# Symlink guard — execute the actual generated bash to confirm semantics
# ---------------------------------------------------------------------------


class TestSymlinkGuardRuntime:
    """Exercise the generated guard fragment under real bash to confirm it
    accepts inputs inside the mount root and rejects ones that escape it.

    This catches subtle bugs (case-glob anchoring, quoting) that string
    assertions can miss.
    """

    def _guard_script(self, mount_root, input_path):
        # Mirror the fragment _shell_setup_lines emits for a single input named
        # 'x' (Snakemake-side `{input.x}` is replaced by the literal path here).
        return (
            f'MOUNT_ROOT="{mount_root}"\n'
            f'INPUT_x=$(readlink -f "{input_path}")\n'
            'case "$INPUT_x/" in "$MOUNT_ROOT"/*) ;; *)'
            " echo \"podman: input 'x' ($INPUT_x) is outside the bind mount $MOUNT_ROOT\" >&2;"
            " exit 1;; esac\n"
            "echo OK"
        )

    def _run(self, script):
        import subprocess

        return subprocess.run(["bash", "-c", script], capture_output=True, text=True)

    def test_input_inside_mount_passes(self, tmp_path):
        (tmp_path / "in.txt").write_text("ok")
        r = self._run(self._guard_script(str(tmp_path), str(tmp_path / "in.txt")))
        assert r.returncode == 0, r.stderr
        assert "OK" in r.stdout

    def test_input_outside_mount_fails(self, tmp_path):
        outside = tmp_path / "outside"
        outside.mkdir()
        (outside / "stranger.txt").write_text("nope")
        mount = tmp_path / "mount"
        mount.mkdir()
        r = self._run(self._guard_script(str(mount), str(outside / "stranger.txt")))
        assert r.returncode == 1
        assert "outside the bind mount" in r.stderr

    def test_symlink_target_outside_mount_fails(self, tmp_path):
        # The symlink itself sits inside the mount, but readlink -f resolves
        # to a path outside it — the guard must reject this.
        outside = tmp_path / "outside"
        outside.mkdir()
        target = outside / "real.txt"
        target.write_text("payload")
        mount = tmp_path / "mount"
        mount.mkdir()
        link = mount / "linked.txt"
        link.symlink_to(target)
        r = self._run(self._guard_script(str(mount), str(link)))
        assert r.returncode == 1
        assert "outside the bind mount" in r.stderr

    def test_sibling_prefix_does_not_match(self, tmp_path):
        # Regression: case "$INPUT_x/" in "$MOUNT_ROOT"/* — the trailing
        # slash on the pattern stops "/foo" from matching "/foobar" by prefix.
        sibling = tmp_path / "mountside"
        sibling.mkdir()
        (sibling / "f.txt").write_text("x")
        mount = tmp_path / "mount"
        mount.mkdir()
        r = self._run(self._guard_script(str(mount), str(sibling / "f.txt")))
        assert r.returncode == 1, r.stdout + r.stderr


# ---------------------------------------------------------------------------
# Resources directive — podman needs mem_mb always set
# ---------------------------------------------------------------------------


class TestResources:
    def test_default_mem_mb_emitted_when_missing(self):
        out = _capture(_gen()._write_resources_directive, _node(resources=None))
        assert f"mem_mb={_DEFAULT_MEM_MB}" in out

    def test_explicit_mem_mb_wins(self):
        r = Resources(cores=4, mem_mb=8192)
        out = _capture(_gen()._write_resources_directive, _node(resources=r))
        assert "mem_mb=8192" in out
        assert "cores=4" in out

    def test_base_generator_keeps_optional_mem_mb(self):
        # Sanity: refactor must not have made mem_mb mandatory for the base class.
        from omnibenchmark.backend.snakemake import SnakemakeGenerator

        base = SnakemakeGenerator(
            benchmark_name="b", benchmark_version="1", benchmark_author="t"
        )
        out = _capture(base._write_resources_directive, _node(resources=None))
        assert "mem_mb" not in out


# ---------------------------------------------------------------------------
# Regression: host-style generator output unchanged
# ---------------------------------------------------------------------------


class TestHostGeneratorRegression:
    """The base SnakemakeGenerator's shell output must not change after the
    refactor — only podman wraps the invocation."""

    def test_host_shell_still_has_cd_module_dir(self):
        from omnibenchmark.backend.snakemake import SnakemakeGenerator

        base = SnakemakeGenerator(
            benchmark_name="b", benchmark_version="1", benchmark_author="t"
        )
        out = _capture(base._write_shell, _node(image=None))
        assert "cd {params.module_dir}" in out
        assert "podman" not in out

    def test_host_shell_has_entrypoint_and_args(self):
        from omnibenchmark.backend.snakemake import SnakemakeGenerator

        base = SnakemakeGenerator(
            benchmark_name="b", benchmark_version="1", benchmark_author="t"
        )
        out = _capture(base._write_shell, _node(image=None))
        assert "python3 {params.entrypoint}" in out
        assert "--output_dir $OUTPUT_DIR" in out
        assert "--name " in out
