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
        f.write("#!/usr/bin/env snakemake\n")
        f.write("#" * 80 + "\n")
        f.write("# Explicit Snakefile for OmniBenchmark\n")
        f.write("#" * 80 + "\n")
        f.write("#\n")
        f.write(f"# Generated: {datetime.now()}\n")
        f.write("#\n")
        f.write("# Benchmark Details:\n")
        f.write(f"#   Name: {self.benchmark_name}\n")
        f.write(f"#   Version: {self.benchmark_version}\n")
        f.write(f"#   Author: {self.benchmark_author}\n")
        f.write("#\n")
        f.write("# This Snakefile was generated from resolved benchmark nodes.\n")
        f.write("# All modules have been cloned and entrypoints dereferenced.\n")
        f.write("# Module checkouts are in .modules/{repo_name}/{commit}/\n")
        f.write("#\n")
        f.write("#" * 80 + "\n")
        f.write("\n")

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

    def _write_shell(self, f: TextIO, node: ResolvedNode):
        """Write the shell: block for a regular node (production execution).

        - Streams stdout/stderr through tee to both terminal and {log}
        - Writes parameters.json next to outputs
        - Uses pre-resolved shebang / interpreter info
        """
        f.write("    shell:\n")
        f.write('        """\n')
        f.write("        mkdir -p $(dirname {output[0]}) $(dirname {log})\n")
        f.write("        OUTPUT_DIR=$(cd $(dirname {output[0]}) && pwd)\n")

        if node.inputs:
            for key in node.inputs.keys():
                f.write(
                    f"        INPUT_{key}=$(cd $(dirname {{input.{key}}}) && pwd)/$(basename {{input.{key}}})\n"
                )

        f.write("        LOG_FILE=$(pwd)/{log}\n")
        f.write('        exec > >(tee "$LOG_FILE") 2>&1\n')
        f.write("        echo '=== Rule: {rule} ==='\n")
        f.write("        echo 'Started:' $(date -Iseconds)\n")
        f.write(
            "        if [ -n \"${{CONDA_PREFIX:-}}\" ]; then echo 'Conda env:' $CONDA_PREFIX; fi\n"
        )
        f.write("        echo '---'\n")

        if node.parameters:
            params_json = json.dumps(node.parameters._params)
            params_json_escaped = params_json.replace("{", "{{").replace("}", "}}")
            params_json_escaped = params_json_escaped.replace("'", "'\\''")
            f.write(
                f"        echo '{params_json_escaped}' > $OUTPUT_DIR/parameters.json\n"
            )
            human_name = _make_human_name(node.parameters)
            hash_folder = f".{node.parameters.hash_short()}"
            f.write(f"        ln -sfn {hash_folder} $OUTPUT_DIR/../{human_name}\n")

        cmd_parts = []
        f.write("        cd {params.module_dir}\n")

        if node.module.has_shebang:
            cmd_parts.append("./{params.entrypoint}")
        else:
            # TODO: do not assume python3, only if .py maybe?
            interpreter = node.module.interpreter or "python3"
            cmd_parts.append(f"{interpreter} {{params.entrypoint}}")

        cmd_parts.append("--output_dir $OUTPUT_DIR")

        # --name is always the current module's own ID (spec §3.5 Reserved Parameters)
        cmd_parts.append(f"--name {node.module_id}")

        if node.inputs:
            for key in node.inputs.keys():
                original_name = node.input_name_mapping.get(key, key)
                cmd_parts.append(f"--{original_name} $INPUT_{key}")

        if node.parameters:
            cmd_parts.append("{params.cli_args}")

        for i, part in enumerate(cmd_parts):
            if i < len(cmd_parts) - 1:
                f.write(f"        {part} \\\\\n")
            else:
                f.write(f"        {part}\n")

        f.write('        """\n')

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

        f.write("    shell:\n")
        f.write('        """\n')
        f.write("        mkdir -p $(dirname {output[0]})\n")
        f.write("        MODULE_DIR={params.module_dir}\n")
        f.write("        OUTPUT_DIR=$(cd $(dirname {output[0]}) && pwd)\n")

        if inputs_by_name:
            input_counter = 0
            for input_name, keys in inputs_by_name.items():
                for key in keys:
                    f.write(
                        f"        INPUT_{input_counter}=$(cd $(dirname {{input.{key}}}) && pwd)/$(basename {{input.{key}}})\n"
                    )
                    input_counter += 1

        if node.module.has_shebang:
            f.write("        $MODULE_DIR/{params.entrypoint} \\\\\n")
        else:
            interpreter = node.module.interpreter or "python3"
            f.write(f"        {interpreter} $MODULE_DIR/{{params.entrypoint}} \\\\\n")

        f.write("            --output_dir $OUTPUT_DIR \\\\\n")
        f.write(f"            --name {node.module_id}")

        if inputs_by_name:
            input_counter = 0
            for input_name, keys in inputs_by_name.items():
                f.write(" \\\\\n")
                f.write(f"            --{input_name}")
                for key in keys:
                    f.write(f" $INPUT_{input_counter}")
                    input_counter += 1

        if node.parameters:
            f.write(" \\\\\n")
            f.write("            {params.cli_args}\n")
        else:
            f.write("\n")

        f.write('        """\n')
