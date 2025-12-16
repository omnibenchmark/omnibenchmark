import logging
import os
import sys
import tempfile


from datetime import datetime
from pathlib import Path
from snakemake.cli import args_to_api as snakemake_cli, parse_args
from typing import TextIO, List, Optional

from omnibenchmark.model import SoftwareBackendEnum

from omnibenchmark.benchmark import BenchmarkNode, BenchmarkExecution
from omnibenchmark.workflow.workflow import WorkflowEngine
from omnibenchmark.workflow.snakemake import rules


# Module includes for each Snakefile
INCLUDES = [
    "utils.smk",
    "rule_start_benchmark.smk",
    "rule_node.smk",
    "rule_all.smk",
]


def _print_ob_flags_help(command_type: str = "benchmark"):
    """Print available ob-specific flags to help users identify typos."""
    print("\n" + "=" * 80, file=sys.stderr)
    print("HINT: The error above is from Snakemake's argument parser.", file=sys.stderr)
    print(
        f"If you intended to use an ob-specific flag, run 'ob run {command_type} --help' for more info.",
        file=sys.stderr,
    )
    print("=" * 80 + "\n", file=sys.stderr)


class SnakemakeEngine(WorkflowEngine):
    """Snakemake implementation of the WorkflowEngine interface."""

    def run_workflow(
        self,
        benchmark: BenchmarkExecution,
        cores: int = 1,
        update: bool = True,
        dryrun: bool = False,
        continue_on_error: bool = False,
        keep_module_logs: bool = False,
        backend: SoftwareBackendEnum = SoftwareBackendEnum.host,
        module_path: str = os.environ.get("MODULEPATH", ""),
        work_dir: Path = Path(os.getcwd()),
        local_timeout: Optional[int] = None,
        debug: bool = False,
        **snakemake_kwargs,
    ) -> bool:
        """
        Serializes & runs benchmark workflow using snakemake.

        Args:
            benchmark (BenchmarkExecution): benchmark to run
            cores (int): number of cores to run. Defaults to 1 core.
            update (bool): run workflow for non-existing outputs / changed nodes only. False means force running workflow from scratch. Default: True
            dryrun (bool): validate the workflow with the benchmark without actual execution. Default: False
            continue_on_error (bool): continue with independent jobs if a job fails. Default: False
            keep_module_logs (bool): keep module-specific log files after execution. Default: False
            backend (SoftwareBackendEnum): which software backend to use when running the workflow. Available: `host`, `docker`, `apptainer`, `conda`, `envmodules`. Default: `host`
            module_path (str): The path where the `envmodules` are located. This path will be searched during the workflow run using `envmodules` backend.
            work_dir (str): working directory. Default: current work directory
            local_timeout (int, optional): timeout, in seconds, when executing locally. It will be passed to the subprocess invocation.
            debug (bool): enable debug output. Default: False
            **snakemake_kwargs: keyword arguments to pass to the snakemake engine

        Returns:
        - Status code (bool) of the workflow run.
        """

        # Get output directory from benchmark context
        out_dir = benchmark.context.out_dir

        # Serialize Snakefile for workflow
        snakefile = self.serialize_workflow(
            benchmark, work_dir, write_to_disk=False, local_timeout=local_timeout
        )

        # Prepare the argv list
        argv = self._prepare_argv(
            snakefile=snakefile,
            cores=cores,
            update=update,
            dryrun=dryrun,
            continue_on_error=continue_on_error,
            keep_module_logs=keep_module_logs,
            backend=backend,
            work_dir=work_dir,
            out_dir=out_dir,
            debug=debug,
            **snakemake_kwargs,
        )

        if module_path:
            os.environ["MODULEPATH"] = module_path

        logging.getLogger("snakemake").setLevel(logging.DEBUG)

        try:
            parser, args = parse_args(argv)
            return snakemake_cli(args, parser)
        except SystemExit as e:
            # Catch argument parsing errors and provide helpful context
            if e.code != 0:
                _print_ob_flags_help(command_type="benchmark")
            raise

    def serialize_workflow(
        self,
        benchmark: BenchmarkExecution,
        output_dir: Path = Path(os.getcwd()),
        write_to_disk=True,
        local_timeout: Optional[int] = None,
    ) -> Path:
        """
        Serializes a Snakefile for the benchmark.

        Args:
            benchmark (BenchmarkExecution): benchmark to serialize
            output_dir (str): output directory for the Snakefile
            write_to_disk (bool): if write_to_disk is True, create Snakefile on disk, else create an in-memory temp file
            local_timeout (int, optional): timeout, in seconds, when executing locally. It will be passed to the subprocess invocation.

        Returns:
        - Snakefile path.
        """
        os.makedirs(output_dir, exist_ok=True)

        benchmark_file = benchmark.get_definition_file()
        name = benchmark.get_benchmark_name()
        version = benchmark.get_benchmark_version()
        author = benchmark.get_benchmark_author()

        if write_to_disk:
            snakefile_path = Path(os.path.join(output_dir, "Snakefile"))
        else:
            temp_file = tempfile.NamedTemporaryFile(
                delete=False, mode="w+", suffix=".smk", encoding="utf-8"
            )
            snakefile_path = Path(temp_file.name)

        # Serialize Snakemake file
        with open(snakefile_path, "w", encoding="utf-8") as f:
            self._write_snakefile_header(f, name, version, author)
            self._write_includes(f, INCLUDES)

            # Load benchmark from yaml file
            f.write(f'benchmark = load("{benchmark_file.as_posix()}", config)\n\n')

            # Create capture all rule
            f.write("node_paths = sorted(benchmark.get_output_paths())\n")
            f.write(
                "mc_paths = sorted(benchmark.get_metric_collector_output_paths())\n"
            )
            f.write("all_paths = node_paths + mc_paths\n")
            f.write(
                "create_all_rule(config, all_paths, aggregate_performance=True)\n\n"
            )

            # Create node rules
            f.write("nodes = benchmark.get_nodes()\n")
            f.write("for node in nodes:\n")
            f.write(
                f"    create_node_rule(node, benchmark, config, {local_timeout})\n\n"
            )

            # Create metric collector rules
            f.write("collectors = benchmark.get_metric_collectors()\n")
            f.write("for collector in collectors:\n")
            f.write(
                "    create_metric_collector_rule(benchmark, collector, config, node_paths)\n\n"
            )

        return snakefile_path

    def run_node_workflow(
        self,
        node: BenchmarkNode,
        input_dir: Path,
        dataset: str,
        cores: int = 1,
        update: bool = True,
        dryrun: bool = False,
        continue_on_error: bool = False,
        keep_module_logs: bool = False,
        backend: SoftwareBackendEnum = SoftwareBackendEnum.host,
        module_path: str = os.environ.get("MODULEPATH", ""),
        work_dir: Path = Path(os.getcwd()),
        local_timeout: Optional[int] = None,
        benchmark_file_path: Optional[Path] = None,
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
            continue_on_error (bool): continue with independent jobs if a job fails. Default: False
            keep_module_logs (bool): keep module-specific log files after execution. Default: False
            backend (SoftwareBackendEnum): which software backend to use when running the workflow. Available: `host`, `docker`, `apptainer`, `conda`, `envmodules`. Default: `host`
            module_path (str): The path where the `envmodules` are located. This path will be searched during the workflow run using `envmodules` backend.
            work_dir (str): working directory. Default: current work directory
            local_timeout (int, optional): timeout, in seconds, when executing locally. It will be passed to the subprocess invocation.
            benchmark_file_path (Optional[Path]): path to benchmark file, if None uses node.get_definition_file()
            **snakemake_kwargs: keyword arguments to pass to the snakemake engine

        Returns:
        - Status code of the workflow run.
        """

        os.makedirs(work_dir, exist_ok=True)

        if benchmark_file_path is None:
            raise ValueError(
                "benchmark_file_path must be provided when node.get_definition_file() returns None"
            )

        # Serialize Snakefile for node workflow
        snakefile = self.serialize_node_workflow(
            node,
            benchmark_file_path=benchmark_file_path,
            output_dir=work_dir,
            write_to_disk=True,
            local_timeout=local_timeout,
        )

        backend_dir = benchmark_file_path.parent.absolute()
        backend_env = node.get_environment_path(
            node.get_software_environment(), backend, backend_dir
        )

        # Prepare the argv list
        argv = self._prepare_argv(
            snakefile=snakefile,
            cores=cores,
            update=update,
            dryrun=dryrun,
            continue_on_error=continue_on_error,
            keep_module_logs=keep_module_logs,
            backend=backend,
            work_dir=work_dir,
            out_dir=work_dir,
            input_dir=input_dir,
            dataset=dataset,
            backend_env=backend_env,
            **snakemake_kwargs,
        )

        if module_path:
            os.environ["MODULEPATH"] = module_path

        try:
            parser, args = parse_args(argv)
            success = snakemake_cli(args, parser)
            return success
        except SystemExit as e:
            # Catch argument parsing errors and provide helpful context
            if e.code != 0:
                _print_ob_flags_help(command_type="module")
            raise

    def serialize_node_workflow(
        self,
        node: BenchmarkNode,
        benchmark_file_path: Optional[Path] = None,
        output_dir: Path = Path(os.getcwd()),
        write_to_disk: bool = True,
        local_timeout: Optional[int] = None,
    ) -> Path:
        """
        Serializes a Snakefile for a benchmark node.

        Args:
            node (BenchmarkNode): benchmark node to serialize
            output_dir (str): output directory for the Snakefile
            write_to_disk (bool): if write_to_disk is True, create Snakefile on disk, else create an in-memory temp file
            local_timeout (int, optional): timeout, in seconds, when executing locally. It will be passed to the subprocess invocation.
            benchmark_file_path (Optional[Path]): path to benchmark file, if None uses node.get_definition_file()

        Returns:
        - Snakefile path.
        """
        if benchmark_file_path is None:
            raise ValueError("benchmark_file_path must be provided")

        os.makedirs(output_dir, exist_ok=True)

        name = node.get_benchmark_name()
        version = node.get_benchmark_version()
        author = node.get_benchmark_author()

        if write_to_disk:
            snakefile_path = Path(os.path.join(output_dir, "Snakefile"))
        else:
            temp_file = tempfile.NamedTemporaryFile(
                delete=False, mode="w+", suffix=".smk", encoding="utf-8"
            )
            snakefile_path = Path(temp_file.name)

        # Serialize Snakemake file
        with open(snakefile_path, "w") as f:
            self._write_snakefile_header(f, name, version, author)
            self._write_includes(f, INCLUDES)

            # Load benchmark from yaml file
            f.write(
                f'node = load_node("{benchmark_file_path.as_posix()}", "{node.get_id()}")\n\n'
            )

            # Create capture all rule
            f.write("input_paths = node.get_input_paths(config)\n")
            f.write("output_paths = node.get_output_paths(config)\n")
            f.write("all_paths = input_paths + output_paths\n\n")
            f.write("create_all_rule(config, all_paths)\n\n")

            # Create node rules
            f.write(f"create_standalone_node_rule(node, config, {local_timeout})\n\n")

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
        keep_module_logs: bool,
        continue_on_error: bool,
        backend: SoftwareBackendEnum,
        work_dir: Path,
        out_dir: Path,
        input_dir: Optional[Path] = None,
        dataset: Optional[str] = None,
        debug: bool = False,
        backend_env: Optional[str] = None,
        **snakemake_kwargs,
    ):
        """Prepare arguments to input to the snakemake cli"""
        argv = [
            "--snakefile",
            str(snakefile),
            "--cores",
            str(cores),
            "--directory",
            work_dir.as_posix(),
            "--config",
            f"dataset={dataset}",
            f"keep_module_logs={keep_module_logs}",
            f"keep_going={continue_on_error}",
            f"backend={backend.value}",
            f"backend_env={backend_env}",
            f"out_dir={out_dir.as_posix()}",
        ]

        if input_dir:
            argv.extend([f"input={input_dir.as_posix()}"])

        if debug:
            argv.extend(["--verbose", "--debug"])

        if update:
            argv.append("-F")

        if dryrun:
            argv.append("--dryrun")

        if continue_on_error:
            argv.append("--keep-going")

        if (
            backend == SoftwareBackendEnum.docker
            or backend == SoftwareBackendEnum.apptainer
        ):
            argv.append("--use-singularity")
            # Bind omnibenchmark package and git repos so they're accessible in the container
            import omnibenchmark
            from omnibenchmark.config import get_git_modules_dir

            omnibenchmark_root = Path(omnibenchmark.__file__).parent.parent.resolve()
            git_modules_dir = get_git_modules_dir().resolve()
            # Use comma-separated bind mounts and set PYTHONPATH for singularity
            bind_paths = f"--bind {omnibenchmark_root}:{omnibenchmark_root},{git_modules_dir}:{git_modules_dir} --env PYTHONPATH={omnibenchmark_root}:$PYTHONPATH"
            argv.extend(["--singularity-args", bind_paths])
        elif backend == SoftwareBackendEnum.envmodules:
            argv.append("--use-envmodules")
        elif backend == SoftwareBackendEnum.conda:
            argv.append("--use-conda")

        for key, value in snakemake_kwargs.items():
            if isinstance(value, bool):
                if value:  # Add flag only if True
                    argv.append(f"--{key}")
            elif isinstance(value, list):
                argv.append(f"--{key}")
                argv.extend(str(v) for v in value)
            elif value is not None:
                argv.extend([f"--{key}", str(value)])

        return argv
