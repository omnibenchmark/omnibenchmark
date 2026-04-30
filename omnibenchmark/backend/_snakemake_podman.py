"""
Podman Snakefile generator.

Wraps each module's host-style invocation in ``podman run`` so cgroup-enforced
CPU/memory limits (derived from the rule's ``resources:`` directive) apply
without depending on Snakemake's container plugins.

Mount model
-----------
The benchmark working directory's *physical* path (``pwd -P``) is bind-mounted
into the container at the same path and used as the workdir.  Inputs and
outputs are then resolved with ``readlink -f`` and asserted to live under that
mount root before the container starts.  The guard turns symlinks pointing
outside the bind mount — common with Snakemake storage plugins, shared module
caches (``ln -s /shared/.modules .modules``), and homedir-symlink CWDs — into
a loud rule failure rather than a confusing ENOENT inside the user's module.

Anything outside the mount root is unsupported.  Use the docker backend (which
goes through Snakemake's container plugin) if you need cross-mount paths.
"""

from typing import Optional, TextIO

from omnibenchmark.backend.snakemake import SnakemakeGenerator
from omnibenchmark.model.resolved import ResolvedModule, ResolvedNode


# Containers always need *some* memory limit for cgroup accounting; this
# default applies only when the benchmark spec leaves mem_mb unset.
_DEFAULT_MEM_MB = 2048


class PodmanSnakemakeGenerator(SnakemakeGenerator):
    """Wrap module invocations in ``podman run`` with resource limits."""

    _default_mem_mb: Optional[int] = _DEFAULT_MEM_MB

    def __init__(self, *args, default_image: Optional[str] = None, **kwargs):
        # default_image is used when a module's resolved environment has no
        # docker reference. If neither is set, _emit_invocation raises.
        super().__init__(*args, **kwargs)
        self._default_image = default_image

    def _write_environment_directive(
        self, f: TextIO, resolved_module: ResolvedModule
    ) -> None:
        # Container is launched explicitly in the shell block; suppress
        # Snakemake's container: directive so --use-singularity has no effect.
        return

    def _shell_setup_lines(self, node: ResolvedNode) -> list:
        """Podman setup: resolve the mount root and validate input symlinks.

        Differs from the base setup in three ways:
          1. ``MOUNT_ROOT`` is the *physical* CWD (``pwd -P``); the bind mount
             uses this so a symlinked CWD doesn't escape the namespace.
          2. ``OUTPUT_DIR`` is also resolved with ``pwd -P``.
          3. Each input is resolved with ``readlink -f`` and a ``case`` guard
             aborts the rule if the resolved path escapes ``MOUNT_ROOT``.
        """
        import json

        from omnibenchmark.backend.snakemake import _make_human_name

        lines: list = [
            "MOUNT_ROOT=$(pwd -P)",
            "mkdir -p {params.output_dir} $(dirname {log})",
            "OUTPUT_DIR=$(cd {params.output_dir} && pwd -P)",
        ]

        for key in node.inputs:
            lines += [
                f"INPUT_{key}=$(readlink -f {{input.{key}}})",
                # case-glob is a POSIX prefix test; the trailing slash on the
                # pattern stops "/foo" from matching "/foobar".
                f'case "$INPUT_{key}/" in "$MOUNT_ROOT"/*) ;; *)'
                f" echo \"podman: input '{key}' ($INPUT_{key}) is outside"
                f' the bind mount $MOUNT_ROOT" >&2; exit 1;; esac',
            ]

        lines += [
            "LOG_FILE=$MOUNT_ROOT/{log}",
            'exec > >(tee "$LOG_FILE") 2>&1',
            "echo '=== Rule: {rule} ==='",
            "echo 'Started:' $(date -Iseconds)",
            "echo '---'",
        ]

        if node.parameters:
            params_json = json.dumps(node.parameters._params)
            params_json_escaped = params_json.replace("{", "{{").replace("}", "}}")
            params_json_escaped = params_json_escaped.replace("'", "'\\''")
            lines += [
                f"echo '{params_json_escaped}' > $OUTPUT_DIR/parameters.json",
                f"ln -sfn .{node.parameters.hash_short()}"
                f" $OUTPUT_DIR/../{_make_human_name(node.parameters)}",
            ]

        return lines

    def _image_for(self, node: ResolvedNode) -> str:
        env = node.module.resolved_environment
        ref = env.reference if env else None
        # Strip any docker:// scheme — podman expects a bare image reference.
        if ref and ref.startswith("docker://"):
            ref = ref[len("docker://") :]
        ref = ref or self._default_image
        if not ref:
            raise ValueError(
                f"podman backend: module '{node.module_id}' has no docker image "
                "and no default_image was provided"
            )
        return ref

    def _emit_invocation(self, node: ResolvedNode, entrypoint: str, args: list) -> list:
        image = self._image_for(node)
        # Single-quoted inner command keeps Snakemake/bash escaping simple. The
        # only character that needs care inside single quotes is "'" itself —
        # none of the generated tokens contain literal single quotes today.
        inner = f"cd {{params.module_dir}} && {entrypoint} " + " ".join(args)
        podman = (
            "podman run --rm"
            " --cpus {resources.cores}"
            " --memory {resources.mem_mb}m"
            " --pull=missing"
            " --volume $MOUNT_ROOT:$MOUNT_ROOT"
            " --workdir $MOUNT_ROOT"
            f" {image}"
            f" sh -c '{inner}'"
        )
        return [podman]
