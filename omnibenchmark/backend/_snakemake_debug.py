"""
Debug (dry-run) Snakefile generator.

Separated from snakemake.py to keep the production module focused on
real execution logic.  Import DebugSnakemakeGenerator from here or via
the backend package.
"""

from collections import defaultdict
from typing import TextIO

from omnibenchmark.backend.snakemake import SnakemakeGenerator, _make_human_name
from omnibenchmark.model.resolved import ResolvedNode


class DebugSnakemakeGenerator(SnakemakeGenerator):
    """Dry-run generator: emits echo/touch commands instead of real execution.

    Useful for testing DAG structure and path resolution without running
    any module code.  Drop-in replacement for SnakemakeGenerator.
    """

    def _write_shell(self, f: TextIO, node: ResolvedNode):
        """Write a debug shell: block that echoes the command and touches outputs."""
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

        f.write(f'        echo "    --name {node.module_id} \\"\n')

        if node.inputs:
            for key in node.inputs.keys():
                original_name = node.input_name_mapping.get(key, key)
                f.write(f'        echo "    --{original_name} {{input.{key}}} \\"\n')

        f.write('        echo "    {params.cli_args}"\n')
        f.write('        echo ""\n')
        f.write("        mkdir -p $(dirname {output[0]})\n")

        if node.parameters:
            human_name = _make_human_name(node.parameters)
            hash_folder = f".{node.parameters.hash_short()}"
            f.write(
                f"        ln -sfn {hash_folder} $(dirname {{output[0]}})/../{human_name}\n"
            )

        for i, _ in enumerate(node.outputs):
            f.write(f"        touch {{output[{i}]}}\n")

        f.write('        """\n')

    def _write_gather_shell(self, f: TextIO, node: ResolvedNode):
        """Write a debug shell: block for gather/collector nodes."""
        inputs_by_name: defaultdict = defaultdict(list)
        if node.inputs and node.input_name_mapping:
            for key in sorted(node.inputs.keys()):
                original_name = node.input_name_mapping.get(key, key)
                inputs_by_name[original_name].append(key)

        node_type = "GATHER STAGE" if node.is_gather else "METRIC COLLECTOR"
        f.write("    shell:\n")
        f.write('        """\n')
        f.write('        echo "=" "=" "=" "=" "=" "=" "=" "="\n')
        f.write(f'        echo "{node_type}: {node.module_id}"\n')
        f.write('        echo "=" "=" "=" "=" "=" "=" "=" "="\n')
        f.write('        echo "Module: {params.module_dir}"\n')
        f.write('        echo "Entrypoint: {params.entrypoint}"\n')

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
