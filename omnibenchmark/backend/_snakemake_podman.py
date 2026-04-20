"""
Podman and denet+podman Snakefile generators.

Instead of delegating container execution to Snakemake (--use-singularity),
these generators wrap the module command in a `podman run` invocation directly
in the shell: block.  This lets us pass cgroup-enforced resource limits
(--cpus, --memory) derived from Snakemake's own resources: directive, without
waiting for upstream Snakemake container execution plugin support.

Hierarchy:
  SnakemakeGenerator
    └── PodmanSnakemakeGenerator      (podman run ...)
          └── DenetPodmanSnakemakeGenerator  (denet podman run ...)
"""

import json
from typing import TextIO

from omnibenchmark.backend.snakemake import SnakemakeGenerator, _make_human_name
from omnibenchmark.model.benchmark import APIVersion
from omnibenchmark.model.resolved import ResolvedNode


class PodmanSnakemakeGenerator(SnakemakeGenerator):
    """Wrap each module's shell command in `podman run` with cgroup resource limits.

    The container image is taken from the module's resolved environment reference
    (same field used by docker/apptainer backends).  Snakemake's container: directive
    is intentionally omitted — podman is invoked directly in the shell: block.
    """

    def _write_benchmark_directive(self, f: TextIO, node) -> None:
        # denet + podman stats provide profiling; Snakemake's psutil-based benchmark
        # monitor interferes with container execution via process-tree traversal.
        pass

    def _write_environment_directive(self, f: TextIO, resolved_module) -> None:
        # Podman images are consumed in _write_shell, not via Snakemake directives.
        pass

    def _podman_prefix(self, node: ResolvedNode) -> str:
        """Build the `podman run ...` prefix for a node's shell command."""
        image = ""
        if node.module.resolved_environment:
            image = node.module.resolved_environment.reference

        cores = "{resources.cores}"
        mem = "{resources.mem_mb}"
        # --pull=never avoids registry traffic during capture; pre-pull images beforehand.
        # --name uses the shell variable set in _write_shell so podman stats can find it.
        return (
            f"podman run --rm"
            f" --name $CONTAINER_NAME"
            f" --cpus {cores}"
            f" --memory {mem}m"
            f" --pull=never"
            f" --volume $(pwd):$(pwd)"
            f" --workdir $(pwd)"
            f" {image}"
        )

    def _write_shell(self, f: TextIO, node: ResolvedNode) -> None:
        lines: list = []

        lines += [
            "mkdir -p {params.output_dir} $(dirname {log})",
            "OUTPUT_DIR=$(cd {params.output_dir} && pwd)",
        ]
        for key in node.inputs:
            lines.append(
                f"INPUT_{key}=$(cd $(dirname {{input.{key}}}) && pwd)"
                f"/$(basename {{input.{key}}})"
            )

        lines += [
            "LOG_FILE=$(pwd)/{log}",
            'exec > >(tee "$LOG_FILE") 2>&1',
            "echo '=== Rule: {rule} ==='",
            "echo 'Started:' $(date -Iseconds)",
            "echo '---'",
            # Unique name scoped to this shell PID so parallel rules don't collide.
            "CONTAINER_NAME={rule}_$$",
            # Poll cgroup-level stats via podman stats (250ms interval) into JSONL.
            # --format json gives clean output without Go-template escaping headaches.
            # Runs in background; killed after the container exits.
            'STATS_FILE="$OUTPUT_DIR/podman_stats.jsonl"',
            '> "$STATS_FILE"',  # truncate so re-runs don't accumulate
            # jq adds ts_ms (epoch ms) to each record.
            # {{...}} is the Snakemake escape for a literal {..} in the shell block.
            # set +e inside the subshell so pipefail from podman stats (container
            # not yet up) does not kill the poller before the container starts.
            "(set +e; while true; do"
            ' podman stats --no-stream --format json "$CONTAINER_NAME" 2>/dev/null'
            " | jq -c --argjson ts \"$(date +%s%3N)\" '.[] + {{ts_ms: $ts}}'"
            ' >> "$STATS_FILE" 2>/dev/null;'
            " sleep 0.25;"
            " done) &",
            "STATS_PID=$!",
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

        if node.module.has_shebang:
            entrypoint_cmd = "./{params.entrypoint}"
        else:
            interpreter = node.module.interpreter or "python3"
            entrypoint_cmd = f"{interpreter} {{params.entrypoint}}"

        name_param = node.module_id
        if self.api_version <= APIVersion.V0_4_0 and node.outputs:
            first_output = node.outputs[0]
            parts = first_output.split("/")
            if len(parts) > 1:
                name_param = parts[1]

        args = [f"--output_dir $OUTPUT_DIR --name {name_param}"]
        for key in node.inputs:
            original_name = node.input_name_mapping.get(key, key)
            args.append(f"--{original_name} $INPUT_{key}")
        if node.parameters:
            args.append("{params.cli_args}")

        # Build sh -c as a single line so the outer shell expands $OUTPUT_DIR/$INPUT_*
        # (double-quoted) and every arg is part of one command, not separate sh commands.
        sh_cmd = (
            f'sh -c "cd {{params.module_dir}} && {entrypoint_cmd} '
            + " ".join(args)
            + '"'
        )
        cmd = [
            f"{self._podman_prefix(node)} \\",
            sh_cmd,
        ]

        lines.extend(cmd)
        lines += [
            "PODMAN_EXIT=$?",
            "kill $STATS_PID 2>/dev/null || true; wait $STATS_PID 2>/dev/null || true",
            "(exit $PODMAN_EXIT)",
        ]
        self._write_shell_lines(f, lines)


class DenetPodmanSnakemakeGenerator(PodmanSnakemakeGenerator):
    """Like PodmanSnakemakeGenerator but prefixes each podman invocation with denet.

    denet runs on the host (eBPF requires host kernel access) and wraps the
    entire podman process, capturing per-PID network counters.  Filter the root
    podman PID from denet's JSON output to get clean container-scoped numbers.
    """

    def __init__(self, *args, denet_output_dir: str = ".denet", **kwargs):
        super().__init__(*args, **kwargs)
        self._denet_output_dir = denet_output_dir

    def _podman_prefix(self, node: ResolvedNode) -> str:
        base = super()._podman_prefix(node)
        rule = self._sanitize_rule_name(node.id)
        out = f"{self._denet_output_dir}/{rule}.json"
        # mkdir ensures the output dir exists; compound so denet only runs after.
        return f"mkdir -p {self._denet_output_dir} && denet --json --enable-ebpf --out {out} run -- {base}"
