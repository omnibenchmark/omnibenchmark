"""
Debug (dry-run) Snakefile generator.

Separated from snakemake.py to keep the production module focused on
real execution logic.  Import DebugSnakemakeGenerator from here or via
the backend package.
"""

from collections import defaultdict
from typing import TextIO

from omnibenchmark.backend.snakemake import (
    SnakemakeGenerator,
    _make_human_name,
)
from omnibenchmark.model.resolved import ResolvedNode


class DebugSnakemakeGenerator(SnakemakeGenerator):
    """Dry-run generator: emits echo/touch commands instead of real execution.

    Useful for testing DAG structure and path resolution without running
    any module code.  Drop-in replacement for SnakemakeGenerator.
    """

    def _write_shell(self, f: TextIO, node: ResolvedNode):
        """Write a debug shell: block that echoes the command and touches outputs."""
        lines: list = [
            'echo "=" "=" "=" "=" "=" "=" "=" "="',
            f'echo "RULE: {node.id}"',
            'echo "=" "=" "=" "=" "=" "=" "=" "="',
            'echo "Module: {params.module_dir}"',
            'echo "Entrypoint: {params.entrypoint}"',
        ]

        if node.inputs:
            lines.append('echo "Inputs:"')
            for key in node.inputs:
                lines.append(f'echo "  {key}: {{input.{key}}}"')

        if node.outputs:
            lines.append('echo "Outputs:"')
            for i in range(len(node.outputs)):
                lines.append(f'echo "  {{output[{i}]}}"')

        if node.parameters:
            lines.append('echo "Parameters: {params.cli_args}"')

        lines += [
            'echo ""',
            'echo "Would execute:"',
            'echo "  cd {params.module_dir} && ./{params.entrypoint} \\"',
            'echo "    --output_dir $(dirname {output[0]}) \\"',
            f'echo "    --name {node.module_id} \\"',
        ]
        for key in node.inputs:
            original_name = node.input_name_mapping.get(key, key)
            lines.append(f'echo "    --{original_name} {{input.{key}}} \\"')
        lines += [
            'echo "    {params.cli_args}"',
            'echo ""',
            "mkdir -p $(dirname {output[0]})",
        ]

        if node.parameters:
            lines.append(
                f"ln -sfn .{node.parameters.hash_short()}"
                f" $(dirname {{output[0]}})/../{_make_human_name(node.parameters)}"
            )

        for i in range(len(node.outputs)):
            lines.append(f"touch {{output[{i}]}}")

        self._write_shell_lines(f, lines)

    def _write_gather_shell(self, f: TextIO, node: ResolvedNode):
        """Write a debug shell: block for gather/collector nodes."""
        inputs_by_name: defaultdict = defaultdict(list)
        if node.inputs and node.input_name_mapping:
            for key in sorted(node.inputs.keys()):
                original_name = node.input_name_mapping.get(key, key)
                inputs_by_name[original_name].append(key)

        node_type = "GATHER STAGE" if node.is_gather else "METRIC COLLECTOR"
        lines: list = [
            'echo "=" "=" "=" "=" "=" "=" "=" "="',
            f'echo "{node_type}: {node.module_id}"',
            'echo "=" "=" "=" "=" "=" "=" "=" "="',
            'echo "Module: {params.module_dir}"',
            'echo "Entrypoint: {params.entrypoint}"',
        ]

        if inputs_by_name:
            lines.append('echo "Inputs (by name):"')
            for input_name, keys in inputs_by_name.items():
                lines.append(f'echo "  --{input_name}:"')
                for key in keys:
                    lines.append(f'echo "    {{input.{key}}}"')

        lines.append('echo "Outputs:"')
        for i in range(len(node.outputs)):
            lines.append(f'echo "  {{output[{i}]}}"')

        if node.parameters:
            lines.append('echo "Parameters: {params.cli_args}"')

        lines += [
            'echo ""',
            'echo "Would execute:"',
            'echo "  cd {params.module_dir} && python3 {params.entrypoint} \\"',
            'echo "    --output_dir $(dirname {output[0]}) \\"',
            f'echo "    --name {node.module_id} \\"',
        ]
        for input_name in inputs_by_name:
            lines.append(f'echo "    --{input_name} [files...] \\"')
        if node.parameters:
            lines.append('echo "    {params.cli_args}"')
        lines += [
            'echo ""',
            "mkdir -p $(dirname {output[0]})",
        ]
        for i in range(len(node.outputs)):
            lines.append(f"touch {{output[{i}]}}")

        self._write_shell_lines(f, lines)
