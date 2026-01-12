"""
Explicit Snakefile generation from resolved nodes.

This module generates explicit, human-readable Snakefiles from ResolvedNode
entities. The generated Snakefiles use shell: directives with direct command
invocation, avoiding runtime complexity.

Design:
- Each ResolvedNode becomes one Snakemake rule
- Rules use shell: directive (not script:)
- Module paths are relative (.repos/{commit}/)
- Parameters are passed via CLI args
- No Python imports or omnibenchmark package dependency
"""

import shutil
from pathlib import Path
from typing import List, TextIO
from datetime import datetime

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

            # Collect dataset values from entrypoint nodes (first stage, no inputs)
            dataset_values = []
            for node in nodes:
                if node.is_entrypoint():
                    # Entrypoint node - module_id is the dataset value
                    if node.module_id not in dataset_values:
                        dataset_values.append(node.module_id)

            # Generate all rule
            self._write_all_rule(f, nodes, collectors, dataset_values)

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
        f.write("# Module checkouts are in .repos/{commit}/\n")
        f.write("#\n")
        f.write("#" * 80 + "\n")
        f.write("\n")

    def _write_node_rule(self, f: TextIO, node: ResolvedNode, debug_mode: bool):
        """Write a rule for a resolved node."""
        f.write(f"# Stage: {node.stage_id}, Module: {node.module_id}\n")
        f.write(f"# Repository: {node.module.repository_url}\n")
        f.write(f"# Commit: {node.module.commit}\n")
        f.write(f"rule {self._sanitize_rule_name(node.id)}:\n")

        # Inputs
        if node.inputs:
            f.write("    input:\n")
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

        # Software environment directive
        self._write_environment_directive(f, node)

        # Shell command
        if debug_mode:
            self._write_debug_shell(f, node)
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
        # Use resolved dataset value if available, otherwise use wildcard
        if node.dataset:
            f.write(f'        echo "    --name {node.dataset} \\"\n')
        else:
            f.write('        echo "    --name {wildcards.dataset} \\"\n')

        if node.inputs:
            for key in node.inputs.keys():
                # Use original name if available, otherwise use sanitized name
                original_name = node.input_name_mapping.get(key, key)
                f.write(f'        echo "    --{original_name} {{input.{key}}} \\"\n')

        f.write('        echo "    {params.cli_args}"\n')
        f.write('        echo ""\n')

        # Create empty output files for dry-run
        f.write("        mkdir -p $(dirname {output[0]})\n")
        for i, _ in enumerate(node.outputs):
            f.write(f"        touch {{output[{i}]}}\n")

        f.write('        """\n')

    def _write_exec_shell(self, f: TextIO, node: ResolvedNode):
        """Write actual execution shell command.

        Uses pre-resolved entrypoint execution info from module resolution:
        - If has_shebang: execute entrypoint directly
        - Otherwise: use inferred interpreter

        Note: After cd .repos/HASH, we use ../../ to get back to out/ directory
        """
        f.write("    shell:\n")
        f.write('        """\n')
        # Create output directory before cd (relative to out/)
        f.write("        mkdir -p $(dirname {output[0]}) && cd {params.module_dir}\n")

        # Use pre-resolved execution info instead of runtime checking
        if node.module.has_shebang:
            # Has shebang, execute directly
            f.write("        ./{params.entrypoint} \\\\\n")
            # Use ../../ to get back to out/ directory from .repos/HASH/
            f.write("            --output_dir $(dirname ../../{output[0]}) \\\\\n")
            # Use resolved dataset value if available, otherwise use wildcard
            if node.dataset:
                f.write(f"            --name {node.dataset}")
            else:
                f.write("            --name {wildcards.dataset}")

            if node.inputs:
                f.write(" \\\\\n")
                for i, key in enumerate(node.inputs.keys()):
                    # Use original name if available, otherwise use sanitized name
                    original_name = node.input_name_mapping.get(key, key)
                    # Prefix input paths with ../../ to get back to out/ directory
                    if i < len(node.inputs) - 1 or node.parameters:
                        f.write(
                            f"            --{original_name} ../../{{input.{key}}} \\\\\n"
                        )
                    else:
                        f.write(f"            --{original_name} ../../{{input.{key}}}")

            if node.parameters:
                f.write(" \\\\\n")
                f.write("            {params.cli_args}\n")
            else:
                f.write("\n")
        else:
            # No shebang, use inferred interpreter
            interpreter = node.module.interpreter or "python3"
            f.write(f"        {interpreter} {{params.entrypoint}} \\\\\n")
            # Use ../../ to get back to out/ directory from .repos/HASH/
            f.write("            --output_dir $(dirname ../../{output[0]}) \\\\\n")
            # Use resolved dataset value if available, otherwise use wildcard
            if node.dataset:
                f.write(f"            --name {node.dataset}")
            else:
                f.write("            --name {wildcards.dataset}")

            if node.inputs:
                f.write(" \\\\\n")
                for i, key in enumerate(node.inputs.keys()):
                    # Use original name if available, otherwise use sanitized name
                    original_name = node.input_name_mapping.get(key, key)
                    # Prefix input paths with ../../ to get back to out/ directory
                    if i < len(node.inputs) - 1 or node.parameters:
                        f.write(
                            f"            --{original_name} ../../{{input.{key}}} \\\\\n"
                        )
                    else:
                        f.write(f"            --{original_name} ../../{{input.{key}}}")

            if node.parameters:
                f.write(" \\\\\n")
                f.write("            {params.cli_args}\n")
            else:
                f.write("\n")

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

        # Software environment directive (for collectors too)
        self._write_environment_directive_for_collector(f, collector)

        # Shell command (simplified for collectors)
        if debug_mode:
            f.write("    shell:\n")
            f.write('        """\n')
            f.write(f'        echo "Metric Collector: {collector.id}"\n')
            f.write('        echo "Module: {params.module_dir}"\n')
            f.write('        echo "Entrypoint: {params.entrypoint}"\n')
            f.write("        mkdir -p $(dirname {output[0]})\n")
            for i, _ in enumerate(collector.outputs):
                f.write(f"        touch {{output[{i}]}}\n")
            f.write('        """\n')
        else:
            f.write("    shell:\n")
            f.write('        """\n')
            f.write("        cd {params.module_dir} && \\\n")
            f.write("        python3 {params.entrypoint} \\\n")
            f.write("            --output {output[0]} \\\n")
            f.write("            {input}\n")
            f.write('        """\n')

        f.write("\n")

    def _write_all_rule(
        self,
        f: TextIO,
        nodes: List[ResolvedNode],
        collectors: List[ResolvedMetricCollector],
        dataset_values: List[str],
    ):
        """
        Write the 'all' rule that depends on all outputs.

        Expands {dataset} wildcards in outputs to concrete dataset values.

        Args:
            f: File handle
            nodes: List of resolved nodes
            collectors: List of resolved metric collectors
            dataset_values: List of dataset values (module IDs from entrypoint nodes)
        """
        f.write("# Target rule: executes entire benchmark\n")
        f.write("rule all:\n")
        f.write("    input:\n")

        # All node outputs (expand {dataset} wildcard)
        for node in nodes:
            for output in node.outputs:
                if "{dataset}" in output:
                    # Expand wildcard to all dataset values
                    for dataset in dataset_values:
                        expanded_output = output.replace("{dataset}", dataset)
                        f.write(f'        "{expanded_output}",\n')
                else:
                    # No wildcard, write as-is
                    f.write(f'        "{output}",\n')

        # All collector outputs (expand {dataset} wildcard)
        for collector in collectors:
            for output in collector.outputs:
                if "{dataset}" in output:
                    # Expand wildcard to all dataset values
                    for dataset in dataset_values:
                        expanded_output = output.replace("{dataset}", dataset)
                        f.write(f'        "{expanded_output}",\n')
                else:
                    # No wildcard, write as-is
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
