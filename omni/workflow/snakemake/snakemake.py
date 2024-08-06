import os
from pathlib import Path
from typing import TextIO, List, Optional

from omni.benchmark import Benchmark, BenchmarkNode
from omni.workflow.workflow import WorkflowEngine
from omni.workflow.snakemake import rules
from snakemake.cli import args_to_api as snakemake_cli, parse_args

from datetime import datetime

# Module includes for each Snakefile
INCLUDES = [
    "utils.smk",
    "rule_start_benchmark.smk",
    "rule_node.smk",
    "rule_all.smk",
]


class SnakemakeEngine(WorkflowEngine):
    """Snakemake implementation of the WorkflowEngine interface."""

    def run_workflow(
        self,
        benchmark: Benchmark,
        cores: int = 1,
        update: bool = True,
        dryrun: bool = False,
        work_dir: Path = Path(os.getcwd()),
        **snakemake_kwargs,
    ) -> bool:
        """
        Serializes & runs benchmark workflow using snakemake.

        Args:
            benchmark (Benchmark): benchmark to run
            cores (int): number of cores to run. Defaults to 1 core.
            update (bool): run workflow for non-existing outputs / changed nodes only. False means force running workflow from scratch. Default: True
            dryrun (bool): validate the workflow with the benchmark without actual execution. Default: False
            work_dir (str): working directory. Default: current work directory
            **snakemake_kwargs: keyword arguments to pass to the snakemake engine

        Returns:
        - Status code (bool) of the workflow run.
        """

        # Serialize Snakefile for workflow
        snakefile = self.serialize_workflow(benchmark, work_dir)

        # Prepare the argv list
        argv = self._prepare_argv(
            snakefile, cores, update, dryrun, work_dir, **snakemake_kwargs
        )

        # Execute snakemake script
        parser, args = parse_args(argv)
        success = snakemake_cli(args, parser)

        return success

    def serialize_workflow(
        self, benchmark: Benchmark, output_dir: Path = Path(os.getcwd())
    ) -> Path:
        """
        Serializes a Snakefile for the benchmark.

        Args:
            benchmark (Benchmark): benchmark to serialize
            output_dir (str): output directory for the Snakefile

        Returns:
        - Snakefile path.
        """
        os.makedirs(output_dir, exist_ok=True)

        benchmark_file = benchmark.get_definition_file()
        name = benchmark.get_benchmark_name()
        version = benchmark.get_benchmark_version()
        author = benchmark.get_benchmark_author()

        # Serialize Snakemake file
        snakefile_path = Path(os.path.join(output_dir, "Snakefile"))
        with open(snakefile_path, "w", encoding="utf-8") as f:
            self._write_snakefile_header(f, name, version, author)
            self._write_includes(f, INCLUDES)

            # Load benchmark from yaml file
            f.write(f'benchmark = load("{benchmark_file.as_posix()}")\n\n')

            # Create capture all rule
            f.write("all_paths = sorted(benchmark.get_output_paths())\n")
            f.write("create_all_rule(all_paths)\n\n")

            # Create node rules
            f.write("nodes = benchmark.get_nodes()\n")
            f.write("for node in nodes:\n")
            f.write("    create_node_rule(node, benchmark)\n\n")

        return snakefile_path

    def run_node_workflow(
        self,
        node: BenchmarkNode,
        input_dir: Path,
        dataset: str,
        cores: int = 1,
        update: bool = True,
        dryrun: bool = False,
        work_dir: Path = Path(os.getcwd()),
        **snakemake_kwargs,
    ) -> bool:
        """
        Serializes & runs benchmark node workflow using snakemake.

        Args:
            node (Benchmark): benchmark node to run
            input_dir (str): directory containing the inputs for the benchmark node
            dataset (str): file names corresponding to the dataset
            cores (int): number of cores to run. Defaults to 1 core.
            update (bool): run workflow for non-existing outputs / changed nodes only. False means force running workflow from scratch. Default: True
            dryrun (bool): validate the workflow with the benchmark without actual execution. Default: False
            work_dir (str): working directory. Default: current work directory
            **snakemake_kwargs: keyword arguments to pass to the snakemake engine

        Returns:
        - Status code of the workflow run.
        """

        os.makedirs(work_dir, exist_ok=True)

        # Serialize Snakefile for node workflow
        snakefile = self.serialize_node_workflow(node, work_dir)

        # Prepare the argv list
        argv = self._prepare_argv(
            snakefile,
            cores,
            update,
            dryrun,
            work_dir,
            input_dir,
            dataset,
            **snakemake_kwargs,
        )

        # Execute snakemake script
        parser, args = parse_args(argv)
        success = snakemake_cli(args, parser)

        return success

    def serialize_node_workflow(
        self, node: BenchmarkNode, output_dir: Path = Path(os.getcwd())
    ) -> Path:
        """
        Serializes a Snakefile for a benchmark node.

        Args:
            node (BenchmarkNode): benchmark node to serialize
            output_dir (str): output directory for the Snakefile

        Returns:
        - Snakefile path.
        """
        os.makedirs(output_dir, exist_ok=True)

        benchmark_file = node.get_definition_file()
        name = node.get_benchmark_name()
        version = node.get_benchmark_version()
        author = node.get_benchmark_author()

        # Serialize Snakemake file
        snakefile_path = Path(os.path.join(output_dir, "Snakefile"))
        with open(snakefile_path, "w") as f:
            self._write_snakefile_header(f, name, version, author)
            self._write_includes(f, INCLUDES)

            # Load benchmark from yaml file
            f.write(
                f'node = load_node("{benchmark_file.as_posix()}", "{node.get_id()}")\n\n'
            )

            # Create capture all rule
            f.write("input_paths = node.get_input_paths(config)\n")
            f.write("output_paths = node.get_output_paths(config)\n")
            f.write("all_paths = input_paths + output_paths\n\n")
            f.write("create_all_rule(all_paths)\n\n")

            # Create node rules
            f.write("create_standalone_node_rule(node, config)\n\n")

        return snakefile_path

    @staticmethod
    def _write_snakefile_header(
        f: TextIO, benchmark_name: str, benchmark_version: str, benchmark_author: str
    ):
        """Write header for the generated Snakefile"""

        f.write("#!/usr/bin/env snakemake -s\n")
        f.write("#############################################\n")
        f.write("# Snakefile for Orchestrating YAML-defined OmniBenchmarks\n")
        f.write("#############################################\n")
        f.write("# \n")
        f.write(
            f"# This Snakefile has been automatically generated on {datetime.now()}\n"
        )
        f.write("# \n")
        f.write("# Benchmark Details:\n")
        f.write(f"# - Name: {benchmark_name}\n")
        f.write(f"# - Version: {benchmark_version}\n")
        f.write("# \n")
        f.write("# Author Information:\n")
        f.write(f"# - Contact: {benchmark_author}\n")
        f.write("#\n")
        f.write("#############################################\n")
        f.write("\n")

    @staticmethod
    def _write_includes(f: TextIO, includes: List[str]):
        """Write includes directive for the generated Snakefile"""

        includes_path = Path(rules.__file__).resolve().parent
        for include in includes:
            include_path = includes_path / include
            f.write(f'include: "{include_path.as_posix()}"\n')

        f.write("\n")

    @staticmethod
    def _prepare_argv(
        snakefile: Path,
        cores: int,
        update: bool,
        dryrun: bool,
        work_dir: Path,
        input_dir: Optional[Path] = None,
        dataset: Optional[str] = None,
        **snakemake_kwargs,
    ):
        """Prepare arguments to input to the snakemake cli"""

        argv = [
            "--snakefile",
            str(snakefile),
            "--cores",
            str(cores),
            "--directory",
            str(work_dir),
            "--config",
            f"input={str(input_dir)}",
            f"dataset={dataset}",
        ]

        if update:
            argv.append("-F")

        if dryrun:
            argv.append("--dryrun")

        for key, value in snakemake_kwargs.items():
            if isinstance(value, bool):
                if value:  # Add flag only if True
                    argv.append(f"--{key}")
            else:
                argv.extend([f"--{key}", str(value)])

        return argv
