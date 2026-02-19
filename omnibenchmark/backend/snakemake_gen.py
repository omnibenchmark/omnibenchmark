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
"""

import json
import platform
import shutil
import sys
from pathlib import Path
from typing import List, Optional, TextIO
from datetime import datetime, timezone

from omnibenchmark.model.resolved import ResolvedNode, ResolvedMetricCollector


class SnakemakeGenerator:
    """
    Generate explicit Snakefiles from resolved nodes.

    The generator creates human-readable Snakefiles that can be inspected
    and executed independently of the omnibenchmark package.
    """

    def __init__(
        self, benchmark_name: str, benchmark_version: str, benchmark_author: str
    ):
        """
        Initialize the generator.

        Args:
            benchmark_name: Benchmark name
            benchmark_version: Benchmark version
            benchmark_author: Benchmark author
        """
        self.benchmark_name = benchmark_name
        self.benchmark_version = benchmark_version
        self.benchmark_author = benchmark_author

    def generate_snakefile(
        self,
        nodes: List[ResolvedNode],
        collectors: List[ResolvedMetricCollector],
        output_path: Path,
        debug_mode: bool = True,
    ):
        """
        Generate a complete Snakefile.

        Args:
            nodes: List of resolved nodes
            collectors: List of resolved metric collectors
            output_path: Path to write Snakefile
            debug_mode: If True, use echo commands instead of actual execution
        """
        with open(output_path, "w") as f:
            self._write_header(f)

            # Generate rule for each node
            for node in nodes:
                self._write_node_rule(f, node, debug_mode)

            # Generate rules for metric collectors
            for collector in collectors:
                self._write_collector_rule(f, collector, debug_mode)

            # Generate all rule
            self._write_all_rule(f, nodes, collectors)

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

    def _write_node_rule(self, f: TextIO, node: ResolvedNode, debug_mode: bool):
        """Write a rule for a resolved node."""
        # Check if this is a collector or gather node
        is_collector = node.stage_id.startswith("_collector_")
        is_gather = node.id.startswith("gather_")

        f.write(f"# Stage: {node.stage_id}, Module: {node.module_id}\n")
        if is_collector:
            f.write("# Type: Metric Collector\n")
        elif is_gather:
            f.write("# Type: Gather Stage\n")
        f.write(f"# Repository: {node.module.repository_url}\n")
        f.write(f"# Commit: {node.module.commit}\n")
        f.write(f"rule {self._sanitize_rule_name(node.id)}:\n")

        # Inputs
        if node.inputs:
            f.write("    input:\n")
            if is_collector:
                # For collectors, write inputs as named entries so we can access them
                # by name in the shell command (e.g., {input.input_0})
                for key, path in node.inputs.items():
                    f.write(f'        {key}="{path}",\n')
            else:
                # Regular nodes: named inputs
                for key, path in node.inputs.items():
                    f.write(f'        {key}="{path}",\n')

        # Outputs
        if node.outputs:
            f.write("    output:\n")
            for output in node.outputs:
                f.write(f'        "{output}",\n')

        # Parameters
        f.write("    params:\n")
        f.write(f'        module_dir="{node.module.module_dir}",\n')
        f.write(f'        entrypoint="{node.module.entrypoint}",\n')

        if node.parameters:
            cli_args = " ".join(node.get_parameter_cli_args())
            f.write(f'        cli_args="{cli_args}",\n')
        else:
            f.write('        cli_args="",\n')

        # Benchmark directive (for performance tracking)
        if not is_collector and not is_gather:
            # Generate benchmark file path in the same directory as the output
            # All path variables (like {dataset}) should have been substituted during node resolution
            if node.outputs:
                first_output = node.outputs[0]
                import os

                # Get the directory and add performance file
                benchmark_dir = (
                    os.path.dirname(first_output) if "/" in first_output else "."
                )
                benchmark_file = f"{benchmark_dir}/clustbench_performance.txt"
                f.write("    benchmark:\n")
                f.write(f'        "{benchmark_file}"\n')

        # Log directive (captures stdout/stderr per rule)
        rule_name = self._sanitize_rule_name(node.id)
        f.write("    log:\n")
        f.write(f'        ".logs/{rule_name}.log"\n')

        # Software environment directive
        self._write_environment_directive(f, node)

        # Resource directives
        self._write_resources_directive(f, node)

        # Shell command
        if debug_mode:
            if is_collector or is_gather:
                self._write_debug_collector_shell_v2(f, node)
            else:
                self._write_debug_shell(f, node)
        else:
            if is_collector or is_gather:
                self._write_exec_collector_shell_v2(f, node)
            else:
                self._write_exec_shell(f, node)

        f.write("\n")

    def _write_debug_shell(self, f: TextIO, node: ResolvedNode):
        """Write debug shell command that echoes what would be executed."""
        f.write("    shell:\n")
        f.write('        """\n')
        f.write('        echo "=" "=" "=" "=" "=" "=" "=" "="\n')
        f.write(f'        echo "RULE: {node.id}"\n')
        f.write('        echo "=" "=" "=" "=" "=" "=" "=" "="\n')
        f.write('        echo "Module: {params.module_dir}"\n')
        f.write('        echo "Entrypoint: {params.entrypoint}"\n')

        if node.inputs:
            f.write('        echo "Inputs:"\n')
            for key in node.inputs.keys():
                f.write(f'        echo "  {key}: {{input.{key}}}"\n')

        if node.outputs:
            f.write('        echo "Outputs:"\n')
            for i, _ in enumerate(node.outputs):
                f.write(f'        echo "  {{output[{i}]}}"\n')

        if node.parameters:
            f.write('        echo "Parameters: {params.cli_args}"\n')

        f.write('        echo ""\n')
        f.write('        echo "Would execute:"\n')
        f.write('        echo "  cd {params.module_dir} && ./{params.entrypoint} \\"\n')
        f.write('        echo "    --output_dir $(dirname {output[0]}) \\"\n')
        # Use dataset from template context
        dataset = (
            node.template_context.provides.get("dataset")
            if node.template_context
            else None
        )
        if dataset:
            f.write(f'        echo "    --name {dataset} \\"\n')

        if node.inputs:
            for key in node.inputs.keys():
                # Use original name if available, otherwise use sanitized name
                original_name = node.input_name_mapping.get(key, key)
                f.write(f'        echo "    --{original_name} {{input.{key}}} \\"\n')

        f.write('        echo "    {params.cli_args}"\n')
        f.write('        echo ""\n')

        # Create empty output files for dry-run
        f.write("        mkdir -p $(dirname {output[0]})\n")

        # Create human-readable symlink for the parameter directory
        if node.parameters:
            human_name = self._make_human_name(node.parameters)
            hash_folder = f".{node.parameters.hash_short()}"
            f.write(
                f"        ln -sfn {hash_folder} $(dirname {{output[0]}})/../{human_name}\n"
            )

        for i, _ in enumerate(node.outputs):
            f.write(f"        touch {{output[{i}]}}\n")

        f.write('        """\n')

    def _write_exec_shell(self, f: TextIO, node: ResolvedNode):
        """Write actual execution shell command.

        Uses pre-resolved entrypoint execution info from module resolution:
        - If has_shebang: execute entrypoint directly
        - Otherwise: use inferred interpreter

        Note: After cd .modules/repo/HASH, we use ../../../ to get back to out/ directory

        Output is tee'd to both stdout (for interactive log viewer) and {log} (for telemetry).
        """
        f.write("    shell:\n")
        f.write('        """\n')
        # Create output directory and log directory before cd (relative to out/)
        f.write("        mkdir -p $(dirname {output[0]}) $(dirname {log})\n")

        # Create human-readable symlink for the parameter directory
        if node.parameters:
            human_name = self._make_human_name(node.parameters)
            hash_folder = f".{node.parameters.hash_short()}"
            f.write(
                f"        ln -sfn {hash_folder} $(dirname {{output[0]}})/../{human_name}\n"
            )

        # Start tee to capture all output to both stdout and per-rule log file
        # Use absolute path for log since we'll cd later
        f.write("        LOG_FILE=$(pwd)/{log}\n")
        f.write('        exec > >(tee "$LOG_FILE") 2>&1\n')

        # Log rule context info
        f.write("        echo '=== Rule: {rule} ==='\n")
        f.write("        echo 'Started:' $(date -Iseconds)\n")
        f.write(
            "        if [ -n \"${{CONDA_PREFIX:-}}\" ]; then echo 'Conda env:' $CONDA_PREFIX; fi\n"
        )
        f.write("        echo '---'\n")

        # Write parameters.json if node has parameters
        if node.parameters:
            # Get JSON string from params (use _params internal OrderedDict)
            import json

            params_json = json.dumps(node.parameters._params)
            # Escape for Snakemake by doubling braces
            params_json_escaped = params_json.replace("{", "{{").replace("}", "}}")
            # Use single quotes for echo to avoid quote escaping issues
            # Need to escape single quotes in the JSON (replace ' with '\\'')
            params_json_escaped = params_json_escaped.replace("'", "'\\''")
            f.write(
                f"        echo '{params_json_escaped}' > $(dirname {{output[0]}})/parameters.json\n"
            )

        # Build the command arguments
        cmd_parts = []

        f.write("        cd {params.module_dir}\n")

        # Use pre-resolved execution info instead of runtime checking
        if node.module.has_shebang:
            # Has shebang, execute directly
            cmd_parts.append("./{params.entrypoint}")
        else:
            # No shebang, use inferred interpreter
            interpreter = node.module.interpreter or "python3"
            cmd_parts.append(f"{interpreter} {{params.entrypoint}}")

        # Use ../../../ to get back to out/ directory from .modules/repo/HASH/
        cmd_parts.append("--output_dir $(dirname ../../../{output[0]})")

        # Use dataset from template context
        dataset = (
            node.template_context.provides.get("dataset")
            if node.template_context
            else None
        )
        if dataset:
            cmd_parts.append(f"--name {dataset}")

        if node.inputs:
            for key in node.inputs.keys():
                # Use original name if available, otherwise use sanitized name
                original_name = node.input_name_mapping.get(key, key)
                # Prefix input paths with ../../../ to get back to out/ directory
                cmd_parts.append(f"--{original_name} ../../../{{input.{key}}}")

        if node.parameters:
            cmd_parts.append("{params.cli_args}")

        # Write command (stdout/stderr already redirected via exec above)
        for i, part in enumerate(cmd_parts):
            if i < len(cmd_parts) - 1:
                f.write(f"        {part} \\\\\n")
            else:
                f.write(f"        {part}\n")

        f.write('        """\n')

    def _write_collector_rule(
        self, f: TextIO, collector: ResolvedMetricCollector, debug_mode: bool
    ):
        """Write a rule for a metric collector."""
        f.write(f"# Metric Collector: {collector.id}\n")
        f.write(f"# Repository: {collector.module.repository_url}\n")
        f.write(f"# Commit: {collector.module.commit}\n")
        f.write(f"rule {self._sanitize_rule_name(collector.id)}:\n")

        # Input patterns
        if collector.input_patterns:
            f.write("    input:\n")
            for pattern in collector.input_patterns:
                f.write(f'        "{pattern}",\n')

        # Outputs
        if collector.outputs:
            f.write("    output:\n")
            for output in collector.outputs:
                f.write(f'        "{output}",\n')

        # Parameters
        f.write("    params:\n")
        f.write(f'        module_dir="{collector.module.module_dir}",\n')
        f.write(f'        entrypoint="{collector.module.entrypoint}",\n')

        if collector.parameters:
            cli_args = " ".join(collector.get_parameter_cli_args())
            f.write(f'        cli_args="{cli_args}",\n')
        else:
            f.write('        cli_args="",\n')

        # Software environment directive (for collectors too)
        self._write_environment_directive_for_collector(f, collector)

        # Shell command
        if debug_mode:
            self._write_debug_collector_shell(f, collector)
        else:
            self._write_exec_collector_shell(f, collector)

        f.write("\n")

    def _write_debug_collector_shell(
        self, f: TextIO, collector: ResolvedMetricCollector
    ):
        """Write debug shell command for a collector."""
        f.write("    shell:\n")
        f.write('        """\n')
        f.write('        echo "=" "=" "=" "=" "=" "=" "=" "="\n')
        f.write(f'        echo "METRIC COLLECTOR: {collector.id}"\n')
        f.write('        echo "=" "=" "=" "=" "=" "=" "=" "="\n')
        f.write('        echo "Module: {params.module_dir}"\n')
        f.write('        echo "Entrypoint: {params.entrypoint}"\n')
        f.write('        echo "Inputs:"\n')
        f.write("        for input_file in {input}; do\n")
        f.write('            echo "  $input_file"\n')
        f.write("        done\n")
        f.write('        echo "Outputs:"\n')
        for i, _ in enumerate(collector.outputs):
            f.write(f'        echo "  {{output[{i}]}}"\n')
        if collector.parameters:
            f.write('        echo "Parameters: {params.cli_args}"\n')
        f.write('        echo ""\n')
        f.write('        echo "Would execute:"\n')
        f.write(
            '        echo "  cd {params.module_dir} && python3 {params.entrypoint} \\"\n'
        )
        f.write('        echo "    --output_dir $(dirname {output[0]}) \\"\n')
        f.write('        echo "    --name ' + f"{collector.id}" + ' \\"\n')
        f.write('        echo "    [input files...] \\"\n')
        if collector.parameters:
            f.write('        echo "    {params.cli_args}"\n')
        f.write('        echo ""\n')
        f.write("        mkdir -p $(dirname {output[0]})\n")
        for i, _ in enumerate(collector.outputs):
            f.write(f"        touch {{output[{i}]}}\n")
        f.write('        """\n')

    def _write_exec_collector_shell(
        self, f: TextIO, collector: ResolvedMetricCollector
    ):
        """Write actual execution shell command for a collector."""
        f.write("    shell:\n")
        f.write('        """\n')
        # Create output directory before cd (relative to out/)
        f.write("        mkdir -p $(dirname {output[0]}) && cd {params.module_dir}\n")

        # Use pre-resolved execution info
        if collector.module.has_shebang:
            # Has shebang, execute directly
            f.write("        ./{params.entrypoint} \\\\\n")
            f.write("            --output_dir $(dirname ../../../{output[0]}) \\\\\n")
            f.write(f"            --name {collector.id}")
        else:
            # No shebang, use inferred interpreter
            interpreter = collector.module.interpreter or "python3"
            f.write(f"        {interpreter} {{params.entrypoint}} \\\\\n")
            f.write("            --output_dir $(dirname ../../../{output[0]}) \\\\\n")
            f.write(f"            --name {collector.id}")

        # Add all input files
        if collector.input_patterns:
            f.write(" \\\\\n")
            # Pass all input files as arguments
            # We need to prefix them with ../../../ to get back to out/ directory
            f.write(
                '            $(printf " %s" $(for f in {input}; do echo "../../../$f"; done))'
            )

        # Add parameters if present
        if collector.parameters:
            f.write(" \\\\\n")
            f.write("            {params.cli_args}\n")
        else:
            f.write("\n")

        f.write('        """\n')

    def _write_debug_collector_shell_v2(self, f: TextIO, node: ResolvedNode):
        """Write debug shell command for a collector/gather node (ResolvedNode version)."""
        # Group inputs by their original names before writing
        from collections import defaultdict

        inputs_by_name = defaultdict(list)
        if node.inputs and node.input_name_mapping:
            for key in sorted(node.inputs.keys()):
                original_name = node.input_name_mapping.get(key, key)
                inputs_by_name[original_name].append(key)

        node_type = (
            "GATHER STAGE" if node.id.startswith("gather_") else "METRIC COLLECTOR"
        )
        f.write("    shell:\n")
        f.write('        """\n')
        f.write('        echo "=" "=" "=" "=" "=" "=" "=" "="\n')
        f.write(f'        echo "{node_type}: {node.module_id}"\n')
        f.write('        echo "=" "=" "=" "=" "=" "=" "=" "="\n')
        f.write('        echo "Module: {params.module_dir}"\n')
        f.write('        echo "Entrypoint: {params.entrypoint}"\n')

        # Show inputs grouped by name
        if inputs_by_name:
            f.write('        echo "Inputs (by name):"\n')
            for input_name in inputs_by_name.keys():
                f.write(f'        echo "  --{input_name}:"\n')
                for key in inputs_by_name[input_name]:
                    f.write(f'        echo "    {{input.{key}}}"\n')

        f.write('        echo "Outputs:"\n')
        for i, _ in enumerate(node.outputs):
            f.write(f'        echo "  {{output[{i}]}}"\n')
        if node.parameters:
            f.write('        echo "Parameters: {params.cli_args}"\n')
        f.write('        echo ""\n')
        f.write('        echo "Would execute:"\n')
        f.write(
            '        echo "  cd {params.module_dir} && python3 {params.entrypoint} \\"\n'
        )
        f.write('        echo "    --output_dir $(dirname {output[0]}) \\"\n')
        f.write(f'        echo "    --name {node.module_id} \\"\n')

        if inputs_by_name:
            for input_name in inputs_by_name.keys():
                f.write(f'        echo "    --{input_name} [files...] \\"\n')

        if node.parameters:
            f.write('        echo "    {params.cli_args}"\n')
        f.write('        echo ""\n')
        f.write("        mkdir -p $(dirname {output[0]})\n")
        for i, _ in enumerate(node.outputs):
            f.write(f"        touch {{output[{i}]}}\n")
        f.write('        """\n')

    def _write_exec_collector_shell_v2(self, f: TextIO, node: ResolvedNode):
        """Write actual execution shell command for a collector node (ResolvedNode version)."""
        # Group inputs by their original names before writing
        from collections import defaultdict

        inputs_by_name = defaultdict(list)
        if node.inputs and node.input_name_mapping:
            for key in sorted(node.inputs.keys()):
                original_name = node.input_name_mapping.get(key, key)
                inputs_by_name[original_name].append(key)

        f.write("    shell:\n")
        f.write('        """\n')
        # Create output directory before cd (relative to out/)
        # For collectors, we store the module path but run from the Snakefile directory
        f.write("        mkdir -p $(dirname {output[0]})\n")
        f.write("        MODULE_DIR={params.module_dir}\n")
        # Convert relative output path to absolute for metric collector
        f.write("        OUTPUT_DIR=$(cd $(dirname {output[0]}) && pwd)\n")

        # Convert input paths to absolute since metric collector scripts may call abspath() from wrong cwd
        # Build absolute path variables for each input
        if inputs_by_name:
            input_counter = 0
            for input_name, keys in inputs_by_name.items():
                for key in keys:
                    f.write(
                        f"        INPUT_{input_counter}=$(cd $(dirname {{input.{key}}}) && pwd)/$(basename {{input.{key}}})\n"
                    )
                    input_counter += 1

        # Use pre-resolved execution info
        if node.module.has_shebang:
            # Has shebang, execute directly
            f.write("        $MODULE_DIR/{params.entrypoint} \\\\\n")
            f.write("            --output_dir $OUTPUT_DIR \\\\\n")
            f.write(f"            --name {node.module_id}")
        else:
            # No shebang, use inferred interpreter
            interpreter = node.module.interpreter or "python3"
            f.write(f"        {interpreter} $MODULE_DIR/{{params.entrypoint}} \\\\\n")
            f.write("            --output_dir $OUTPUT_DIR \\\\\n")
            f.write(f"            --name {node.module_id}")

        # Add input files grouped by their original names using the absolute path variables
        if inputs_by_name:
            input_counter = 0
            for input_name, keys in inputs_by_name.items():
                f.write(" \\\\\n")
                f.write(f"            --{input_name}")
                for key in keys:
                    f.write(f" $INPUT_{input_counter}")
                    input_counter += 1

        # Add parameters if present
        if node.parameters:
            f.write(" \\\\\n")
            f.write("            {params.cli_args}\n")
        else:
            f.write("\n")

        f.write('        """\n')

    def _write_all_rule(
        self,
        f: TextIO,
        nodes: List[ResolvedNode],
        collectors: List[ResolvedMetricCollector],
    ):
        """
        Write the 'all' rule that depends on all outputs.

        All output paths are already fully resolved (no wildcards).

        Args:
            f: File handle
            nodes: List of resolved nodes
            collectors: List of resolved metric collectors
        """
        f.write("# Target rule: executes entire benchmark\n")
        f.write("rule all:\n")
        f.write("    input:\n")

        for node in nodes:
            for output in node.outputs:
                f.write(f'        "{output}",\n')

        for collector in collectors:
            for output in collector.outputs:
                f.write(f'        "{output}",\n')

        f.write("    default_target: True\n")
        f.write("\n")

    def _write_environment_directive(self, f: TextIO, node: ResolvedNode):
        """
        Write software environment directive for a node.

        This adds the appropriate Snakemake directive (conda:, container:, envmodules:)
        based on the resolved environment.
        """
        if not node.module.resolved_environment:
            # No environment specified (host backend or missing)
            return

        env = node.module.resolved_environment

        if env.backend_type.value == "conda":
            # Conda environment - use conda: directive
            f.write(f'    conda: "{env.reference}"\n')

        elif env.backend_type.value == "apptainer":
            # Apptainer/Singularity container directive
            # Snakemake uses: container: "path/to/image.sif"
            # or: container: "docker://..." or "oras://..."
            f.write(f'    container: "{env.reference}"\n')

        elif env.backend_type.value == "envmodules":
            # Envmodules directive
            # Snakemake uses: envmodules: "module_name"
            # This requires Snakemake's --use-envmodules flag at execution time
            f.write(f'    envmodules: "{env.reference}"\n')

        elif env.backend_type.value == "docker":
            # Docker container directive
            # Snakemake uses: container: "docker://..."
            f.write(f'    container: "{env.reference}"\n')

        # host backend doesn't need any directive

    def _write_environment_directive_for_collector(
        self, f: TextIO, collector: ResolvedMetricCollector
    ):
        """
        Write software environment directive for a metric collector.

        Similar to _write_environment_directive but for collectors.
        """
        if not collector.module.resolved_environment:
            return

        env = collector.module.resolved_environment

        if env.backend_type.value == "conda":
            f.write(f'    conda: "{env.reference}"\n')
        elif env.backend_type.value == "apptainer":
            f.write(f'    container: "{env.reference}"\n')
        elif env.backend_type.value == "envmodules":
            f.write(f'    envmodules: "{env.reference}"\n')
        elif env.backend_type.value == "docker":
            f.write(f'    container: "{env.reference}"\n')

    def _write_resources_directive(self, f: TextIO, node: ResolvedNode):
        """
        Write resource directives for a node.

        Generates Snakemake resources: directive based on the node's resource requirements.
        Default: 2 cores if no resources specified.

        Hierarchy (most specific wins):
        1. Module-level resources (stored in node.resources)
        2. Stage-level resources (stored in node.resources)
        3. Global default (2 cores)
        """
        # Default resources
        default_cores = 2

        if node.resources is None:
            # No resources specified, use default
            f.write("    resources:\n")
            f.write(f"        cores={default_cores}\n")
            return

        # Resources specified - write them out
        f.write("    resources:\n")

        # Write cores (use default if not specified)
        cores = getattr(node.resources, "cores", None) or default_cores
        f.write(f"        cores={cores}")

        # Write optional resources if specified
        if hasattr(node.resources, "mem_mb") and node.resources.mem_mb:
            f.write(f",\n        mem_mb={node.resources.mem_mb}")

        if hasattr(node.resources, "disk_mb") and node.resources.disk_mb:
            f.write(f",\n        disk_mb={node.resources.disk_mb}")

        if hasattr(node.resources, "runtime") and node.resources.runtime:
            f.write(f",\n        runtime={node.resources.runtime}")

        f.write("\n")

    def _make_human_name(self, params) -> str:
        """Create a human-readable name from parameters for symlinks.

        Replicates SymlinkManager._make_human_name logic so symlinks can be
        emitted directly in generated shell commands.
        """
        parts = [f"{k}-{v}" for k, v in params.items()]
        name = "_".join(parts)
        for unsafe in ["/", "\\", ":", "*", "?", '"', "<", ">", "|", " "]:
            name = name.replace(unsafe, "_")
        if len(name) > 255:
            name = name[:246] + "_" + params.hash_short()
        return name

    def _sanitize_rule_name(self, node_id: str) -> str:
        """
        Sanitize node ID to valid Snakemake rule name.

        Snakemake rule names must be valid Python identifiers.
        """
        # Replace hyphens and dots with underscores
        name = node_id.replace("-", "_").replace(".", "_")

        # Ensure it starts with a letter
        if name and not name[0].isalpha():
            name = "rule_" + name

        return name


def save_metadata(
    benchmark_yaml_path: Path,
    output_dir: Path,
    nodes: List[ResolvedNode],
    collectors: List[ResolvedMetricCollector],
):
    """
    Save benchmark metadata to output directory.

    Creates:
    - out/.metadata/benchmark.yaml (copy of original)
    - out/.metadata/modules.txt (list of modules with URLs and commits)

    Args:
        benchmark_yaml_path: Path to original benchmark YAML
        output_dir: Output directory
        nodes: Resolved nodes
        collectors: Resolved metric collectors
    """
    metadata_dir = output_dir / ".metadata"
    metadata_dir.mkdir(parents=True, exist_ok=True)

    # Copy original benchmark YAML
    shutil.copy(benchmark_yaml_path, metadata_dir / "benchmark.yaml")

    # Write modules list
    with open(metadata_dir / "modules.txt", "w") as f:
        f.write("# Modules used in this benchmark\n")
        f.write(f"# Generated: {datetime.now()}\n")
        f.write("#\n")
        f.write("# Format: stage/module - repository@commit\n")
        f.write("#\n\n")

        # Collect unique modules
        seen_modules = set()

        for node in nodes:
            module_key = (node.module.repository_url, node.module.commit)
            if module_key not in seen_modules:
                seen_modules.add(module_key)
                f.write(f"{node.stage_id}/{node.module_id}:\n")
                f.write(f"  Repository: {node.module.repository_url}\n")
                f.write(f"  Commit: {node.module.commit}\n")
                f.write(f"  Module dir: {node.module.module_dir}\n")
                f.write(f"  Entrypoint: {node.module.entrypoint}\n")
                f.write("\n")

        for collector in collectors:
            module_key = (collector.module.repository_url, collector.module.commit)
            if module_key not in seen_modules:
                seen_modules.add(module_key)
                f.write(f"metric_collector/{collector.id}:\n")
                f.write(f"  Repository: {collector.module.repository_url}\n")
                f.write(f"  Commit: {collector.module.commit}\n")
                f.write(f"  Module dir: {collector.module.module_dir}\n")
                f.write(f"  Entrypoint: {collector.module.entrypoint}\n")
                f.write("\n")


def write_run_manifest(
    output_dir: Path,
    run_id: Optional[str] = None,
) -> None:
    """Write out/.metadata/manifest.json with a stable run UUID and host metadata.

    The *run_id* is the same UUID that the telemetry emitter uses as its
    OTLP trace_id, so every piece of provenance (telemetry.jsonl, manifest.json,
    log files) can be correlated by that single identifier.  When telemetry is
    disabled a fresh UUID is generated here.

    Fields written:
      run_id        – UUID4 string (from telemetry trace_id when available)
      timestamp     – ISO-8601 UTC timestamp of when the manifest was written
      hostname      – machine hostname
      platform      – OS platform string (e.g. "linux")
      os            – full OS release string (e.g. "Linux 6.x #1 SMP …")
      kernel        – kernel release (uname -r equivalent)
      cpu_count     – logical CPU count visible to the process
      cpu_model     – CPU model string where available (Linux /proc/cpuinfo)
      memory_total_mb  – total physical RAM in MiB (where available)
      python_version   – Python version string (e.g. "3.11.8")
      python_executable – path to the Python interpreter used to run the program
      gpu_devices      – list of NVIDIA GPU dicts {index, name, memory_total_mb}
                         from nvidia-smi; null if nvidia-smi is absent or fails
    """
    import uuid

    metadata_dir = output_dir / ".metadata"
    metadata_dir.mkdir(parents=True, exist_ok=True)

    if run_id is None:
        run_id = str(uuid.uuid4())

    uname = platform.uname()

    cpu_count = None
    try:
        import os

        cpu_count = os.cpu_count()
    except Exception:
        pass

    cpu_model = None
    try:
        if sys.platform.startswith("linux"):
            with open("/proc/cpuinfo") as fh:
                for line in fh:
                    if line.startswith("model name"):
                        cpu_model = line.split(":", 1)[1].strip()
                        break
        elif sys.platform == "darwin":
            import subprocess

            result = subprocess.run(
                ["sysctl", "-n", "machdep.cpu.brand_string"],
                capture_output=True,
                text=True,
                timeout=2,
            )
            if result.returncode == 0:
                cpu_model = result.stdout.strip() or None
    except Exception:
        pass

    memory_total_mb = None
    try:
        if sys.platform.startswith("linux"):
            with open("/proc/meminfo") as fh:
                for line in fh:
                    if line.startswith("MemTotal:"):
                        kb = int(line.split()[1])
                        memory_total_mb = kb // 1024
                        break
        elif sys.platform == "darwin":
            import subprocess

            result = subprocess.run(
                ["sysctl", "-n", "hw.memsize"],
                capture_output=True,
                text=True,
                timeout=2,
            )
            if result.returncode == 0:
                memory_total_mb = int(result.stdout.strip()) // (1024 * 1024)
    except Exception:
        pass

    # GPU devices via nvidia-smi (NVIDIA) and rocm-smi (AMD); both optional.
    # nvidia-smi --query-gpu=... --format=csv,noheader,nounits is supported on
    # all modern driver versions and needs no extra Python dependencies.
    gpu_devices = None
    try:
        import subprocess as _sp

        result = _sp.run(
            [
                "nvidia-smi",
                "--query-gpu=index,name,memory.total",
                "--format=csv,noheader,nounits",
            ],
            capture_output=True,
            text=True,
            timeout=5,
        )
        if result.returncode == 0 and result.stdout.strip():
            gpu_devices = []
            for line in result.stdout.strip().splitlines():
                parts = [p.strip() for p in line.split(",")]
                entry = {"index": int(parts[0]), "name": parts[1]}
                if len(parts) >= 3:
                    try:
                        entry["memory_total_mb"] = int(parts[2])
                    except ValueError:
                        pass
                gpu_devices.append(entry)
    except Exception:
        pass

    manifest = {
        "run_id": run_id,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "hostname": uname.node,
        "platform": sys.platform,
        "os": f"{uname.system} {uname.release} {uname.version}".strip(),
        "kernel": uname.release,
        "cpu_count": cpu_count,
        "cpu_model": cpu_model,
        "memory_total_mb": memory_total_mb,
        "python_version": platform.python_version(),
        "python_executable": sys.executable,
        "gpu_devices": gpu_devices,
    }

    manifest_path = metadata_dir / "manifest.json"
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)
        fh.write("\n")
