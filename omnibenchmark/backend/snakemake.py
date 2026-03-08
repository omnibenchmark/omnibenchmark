"""
Explicit Snakefile generation from resolved nodes.

This module generates explicit, human-readable Snakefiles from ResolvedNode
entities. The generated Snakefiles use shell: directives with direct command
invocation, avoiding runtime complexity.

Design:
- Each ResolvedNode becomes one Snakemake rule
- Rules use shell: directive (not script:)
- Module paths are relative (.modules/{commit}/)
- Parameters are passed via CLI args
- No Python imports or omnibenchmark package dependency

The production generator is SnakemakeGenerator.  For dry-run/debug mode
see omnibenchmark.backend._snakemake_debug.DebugSnakemakeGenerator.
"""

import json
import os
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import List, TextIO

from omnibenchmark.model.benchmark import APIVersion
from omnibenchmark.model.resolved import ResolvedModule, ResolvedNode

_UNSAFE_CHARS = ["/", "\\", ":", "*", "?", '"', "<", ">", "|", " "]

# Indentation constants for generated Snakemake syntax.
# Shell block content lives at 8 spaces (2 levels × 4 spaces).
_RULE_INDENT = "    "
_SHELL_INDENT = "        "


def _bash_var(name: str) -> str:
    """Turn an input flag name into a valid bash array variable name.

    >>> _bash_var("data.raw")
    '_data_raw'
    """
    return "_" + name.replace(".", "_").replace("-", "_")


def _make_human_name(params) -> str:
    """Create a human-readable symlink name from a Params object.

    Joins key-value pairs with underscores and strips filesystem-unsafe
    characters.  Falls back to appending a short hash if the result exceeds
    255 characters (the typical filesystem limit).
    """
    parts = [f"{k}-{v}" for k, v in params.items()]
    name = "_".join(parts)
    for char in _UNSAFE_CHARS:
        name = name.replace(char, "_")
    if len(name) > 255:
        name = name[:246] + "_" + params.hash_short()
    return name


class SnakemakeGenerator:
    """Generate explicit Snakefiles from resolved nodes (production mode).

    Each node becomes one Snakemake rule whose shell: block executes the
    module entrypoint directly.  Override _write_shell / _write_gather_shell
    to change the generated shell commands (see DebugSnakemakeGenerator).
    """

    def __init__(
        self,
        benchmark_name: str,
        benchmark_version: str,
        benchmark_author: str,
        api_version: APIVersion = APIVersion.V0_4_0,
    ):
        self.benchmark_name = benchmark_name
        self.benchmark_version = benchmark_version
        self.benchmark_author = benchmark_author
        self.api_version = api_version

    def generate_snakefile(self, nodes: List[ResolvedNode], output_path: Path):
        """Write a complete Snakefile to *output_path*."""
        with open(output_path, "w") as f:
            self._write_header(f)
            for node in nodes:
                self._write_node_rule(f, node)
            self._write_all_rule(f, nodes)

    # ------------------------------------------------------------------
    # Rule structure (shared between production and debug)
    # ------------------------------------------------------------------

    def _write_header(self, f: TextIO):
        """Write Snakefile header."""
        lines = [
            "#!/usr/bin/env snakemake",
            f"# OmniBenchmark: {self.benchmark_name} v{self.benchmark_version} ({self.benchmark_author})",
            f"# Generated: {datetime.now()}",
            "# Auto-generated — do not edit by hand.",
            "",
        ]
        f.write("\n".join(lines))

    def _write_node_rule(self, f: TextIO, node: ResolvedNode):
        """Write a complete rule block for one resolved node."""
        is_collector = node.is_collector
        is_gather = node.is_gather

        f.write(f"# Stage: {node.stage_id}, Module: {node.module_id}\n")
        if is_collector:
            f.write("# Type: Metric Collector\n")
        elif is_gather:
            f.write("# Type: Gather Stage\n")
        f.write(f"# Repository: {node.module.repository_url}\n")
        f.write(f"# Commit: {node.module.commit}\n")
        f.write(f"rule {self._sanitize_rule_name(node.id)}:\n")

        if node.inputs:
            f.write("    input:\n")
            for key, path in node.inputs.items():
                f.write(f'        {key}="{path}",\n')

        if node.outputs:
            f.write("    output:\n")
            for output in node.outputs:
                f.write(f'        "{output}",\n')

        f.write("    params:\n")
        f.write(f'        module_dir="{node.module.module_dir}",\n')
        f.write(f'        entrypoint="{node.module.entrypoint}",\n')
        if node.outputs:
            output_dir = os.path.dirname(node.outputs[0])
            f.write(f'        output_dir="{output_dir}",\n')
        if node.parameters:
            cli_args = " ".join(node.get_parameter_cli_args())
            f.write(f'        cli_args="{cli_args}",\n')
        else:
            f.write('        cli_args="",\n')

        # benchmark: directive (performance tracking, skipped for aggregate nodes)
        if not is_collector and not is_gather and node.outputs:
            first_output = node.outputs[0]
            benchmark_dir = (
                os.path.dirname(first_output) if "/" in first_output else "."
            )
            if self.api_version >= APIVersion.V0_5_0:
                benchmark_file = f"{benchmark_dir}/performance.txt"
            else:
                # COMPAT(0.4): collector R scripts scan for
                # "{dataset}_performance.txt".  parts[1] is the root module name.
                # Remove when 0.4 compat is dropped.
                parts = first_output.split("/")
                dataset_name = parts[1] if len(parts) > 2 else parts[0].split(".")[0]
                benchmark_file = f"{benchmark_dir}/{dataset_name}_performance.txt"
            f.write("    benchmark:\n")
            f.write(f'        "{benchmark_file}"\n')

        rule_name = self._sanitize_rule_name(node.id)
        f.write("    log:\n")
        f.write(f'        ".logs/{rule_name}.log"\n')

        self._write_environment_directive(f, node.module)
        self._write_resources_directive(f, node)

        if is_collector or is_gather:
            self._write_gather_shell(f, node)
        else:
            self._write_shell(f, node)

        f.write("\n")

    def _write_all_rule(self, f: TextIO, nodes: List[ResolvedNode]):
        """Write the 'all' rule that depends on all node outputs."""
        f.write("# Target rule: executes entire benchmark\n")
        f.write("rule all:\n")
        f.write("    input:\n")
        for node in nodes:
            for output in node.outputs:
                f.write(f'        "{output}",\n')
        f.write("    default_target: True\n")
        f.write("\n")

    def _write_environment_directive(self, f: TextIO, resolved_module: ResolvedModule):
        """Write the appropriate Snakemake environment directive for a resolved module.

        Adds conda:, container:, or envmodules: based on the resolved environment.
        No-op for host backend or when no environment is specified.
        """
        if not resolved_module.resolved_environment:
            return

        env = resolved_module.resolved_environment
        if env.backend_type.value == "conda":
            f.write(f'    conda: "{env.reference}"\n')
        elif env.backend_type.value in ("apptainer", "docker"):
            f.write(f'    container: "{env.reference}"\n')
        elif env.backend_type.value == "envmodules":
            f.write(f'    envmodules: "{env.reference}"\n')
        # host backend needs no directive

    def _write_resources_directive(self, f: TextIO, node: ResolvedNode):
        """Write Snakemake resources: directive for a node.

        Falls back to 2 cores when no resources are specified.
        """
        default_cores = 2
        f.write("    resources:\n")

        r = node.resources
        cores = (r.cores if r is not None else None) or default_cores
        f.write(f"        cores={cores}")

        if r is not None:
            if r.mem_mb:
                f.write(f",\n        mem_mb={r.mem_mb}")
            if r.disk_mb:
                f.write(f",\n        disk_mb={r.disk_mb}")
            if r.runtime:
                f.write(f",\n        runtime={r.runtime}")
            if r.gpu:
                f.write(f",\n        nvidia_gpu={r.gpu}")

        f.write("\n")

    def _sanitize_rule_name(self, node_id: str) -> str:
        """Sanitize node ID to a valid Snakemake rule name (Python identifier)."""
        name = node_id.replace("-", "_").replace(".", "_")
        if name and not name[0].isalpha():
            name = "rule_" + name
        return name

    # ------------------------------------------------------------------
    # Shell writers — override in subclasses to change generated commands
    # ------------------------------------------------------------------

    def _write_shell_lines(self, f: TextIO, lines: list) -> None:
        """Write a shell: block from a list of content lines.

        Handles the directive header, triple-quote delimiters, and uniform
        8-space indentation so callers only need to think about content.
        """
        f.write(f"{_RULE_INDENT}shell:\n")
        f.write(f'{_SHELL_INDENT}"""\n')
        for line in lines:
            f.write(f"{_SHELL_INDENT}{line}\n")
        f.write(f'{_SHELL_INDENT}"""\n')

    def _write_shell(self, f: TextIO, node: ResolvedNode):
        """Write the shell: block for a regular node (production execution).

        - Streams stdout/stderr through tee to both terminal and {log}
        - Writes parameters.json next to outputs
        - Uses pre-resolved shebang / interpreter info
        """
        lines: list = []

        # Resolve output and log directories to absolute paths before cd-ing
        # into the module directory.
        lines += [
            "mkdir -p {params.output_dir} $(dirname {log})",
            "OUTPUT_DIR=$(cd {params.output_dir} && pwd)",
        ]
        for key in node.inputs:
            lines.append(
                f"INPUT_{key}=$(cd $(dirname {{input.{key}}}) && pwd)"
                f"/$(basename {{input.{key}}})"
            )

        # Redirect stdout/stderr through tee so both the terminal and the log
        # file receive all output.
        lines += [
            "LOG_FILE=$(pwd)/{log}",
            'exec > >(tee "$LOG_FILE") 2>&1',
            "echo '=== Rule: {rule} ==='",
            "echo 'Started:' $(date -Iseconds)",
            "echo '---'",
        ]

        # Write parameters.json and a human-readable symlink to the hash folder.
        if node.parameters:
            params_json = json.dumps(node.parameters._params)
            params_json_escaped = params_json.replace("{", "{{").replace("}", "}}")
            params_json_escaped = params_json_escaped.replace("'", "'\\''")
            lines += [
                f"echo '{params_json_escaped}' > $OUTPUT_DIR/parameters.json",
                f"ln -sfn .{node.parameters.hash_short()}"
                f" $OUTPUT_DIR/../{_make_human_name(node.parameters)}",
            ]

        lines.append("cd {params.module_dir}")

        # Build command with backslash-continuations between parts.
        # --name is always the current module's own ID (spec §3.5 Reserved Parameters)
        # TODO: do not assume python3, only if .py maybe?
        if node.module.has_shebang:
            cmd = ["./{params.entrypoint}"]
        else:
            interpreter = node.module.interpreter or "python3"
            cmd = [f"{interpreter} {{params.entrypoint}}"]

        # Determine the value for --name parameter
        # BUG WORKAROUND (API <= 0.4): Legacy modules have a design flaw where they
        # use the --name parameter to construct output filenames ({name}_data.json),
        # but the benchmark spec expects outputs to be named after datasets
        # ({dataset}_data.json). This creates a mismatch when module_id != dataset_id.
        #
        # For API <= 0.4 ONLY, we work around this by passing the dataset name
        # (extracted from the output path structure) as --name instead of module_id.
        # This maintains backward compatibility with existing modules.
        #
        # For API >= 0.5, modules MUST use output_dir correctly and not depend on
        # --name for filename construction. Remove this workaround when 0.4 support ends.
        name_param = node.module_id
        if self.api_version <= APIVersion.V0_4_0 and node.outputs:
            first_output = node.outputs[0]
            parts = first_output.split("/")
            # Output path structure: data/{dataset}/.../[methods/{module}/...]/filename
            # We extract parts[1] which is the dataset identifier
            # Examples:
            #   data/D1/.94686c86/D1_data.json -> D1
            #   data/D1/.94686c86/methods/M1/.68f00a0d/D1_data.json -> D1
            if len(parts) > 1:
                dataset_name = parts[1]
                name_param = dataset_name

        cmd += ["--output_dir $OUTPUT_DIR", f"--name {name_param}"]
        for key in node.inputs:
            original_name = node.input_name_mapping.get(key, key)
            cmd.append(f"--{original_name} $INPUT_{key}")
        if node.parameters:
            cmd.append("{params.cli_args}")

        for part in cmd[:-1]:
            lines.append(f"{part} \\\\")
        lines.append(cmd[-1])

        self._write_shell_lines(f, lines)

    def _write_gather_shell(self, f: TextIO, node: ResolvedNode):
        """Write the shell: block for a gather/collector node (production execution).

        Uses absolute paths for inputs and outputs so that collector scripts
        can safely call os.path.abspath() from any working directory.
        """
        inputs_by_name: defaultdict = defaultdict(list)
        if node.inputs and node.input_name_mapping:
            for key in sorted(node.inputs.keys()):
                original_name = node.input_name_mapping.get(key, key)
                inputs_by_name[original_name].append(key)

        lines: list = [
            "mkdir -p {params.output_dir}",
            "MODULE_DIR={params.module_dir}",
            "OUTPUT_DIR=$(cd {params.output_dir} && pwd)",
        ]

        # Resolve each input group to absolute paths into a bash array.
        # One loop per distinct flag name — no matter how many files per flag.
        for input_name, keys in inputs_by_name.items():
            var = _bash_var(input_name)
            refs = " ".join(f"{{input.{key}}}" for key in keys)
            lines += [
                f"{var}=()",
                f"for _f in {refs}; do",
                f'    {var}+=("$(cd $(dirname $_f) && pwd)/$(basename $_f)")',
                "done",
            ]

        # Build the multi-line command.  Continuation args are indented 4 extra
        # spaces relative to the entrypoint line so they visually align.
        if node.module.has_shebang:
            lines.append("$MODULE_DIR/{params.entrypoint} \\\\")
        else:
            interpreter = node.module.interpreter or "python3"
            lines.append(f"{interpreter} $MODULE_DIR/{{params.entrypoint}} \\\\")

        # Gather args: each on its own continuation line, 4 extra spaces.
        # Array expansion: ${{var[@]}} in the Snakefile so Snakemake unescapes {{ → {
        # and }} → } leaving bash with ${var[@]}.
        gather_args = ["--output_dir $OUTPUT_DIR", f"--name {node.module_id}"]
        for input_name in inputs_by_name:
            var = _bash_var(input_name)
            array_ref = '"${{' + var + '[@]}}"'
            gather_args.append(f"--{input_name} {array_ref}")
        if node.parameters:
            gather_args.append("{params.cli_args}")

        for arg in gather_args[:-1]:
            lines.append(f"    {arg} \\\\")
        lines.append(f"    {gather_args[-1]}")

        self._write_shell_lines(f, lines)
