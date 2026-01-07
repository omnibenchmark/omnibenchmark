"""cli commands related to benchmark/module execution and start"""

import os
import shutil
import sys
from itertools import chain
from pathlib import Path

import click
import humanfriendly

from omnibenchmark.benchmark import BenchmarkExecution
from omnibenchmark.cli.utils.args import parse_extra_args
from omnibenchmark.cli.utils.logging import logger
from omnibenchmark.cli.error_formatting import pretty_print_parse_error
from omnibenchmark.model.validation import BenchmarkParseError
from omnibenchmark.remote.storage import get_storage_from_benchmark
from omnibenchmark.remote.storage import remote_storage_snakemake_args
from omnibenchmark.workflow.snakemake import SnakemakeEngine
from omnibenchmark.workflow.workflow import WorkflowEngine


@click.command(
    name="run",
    context_settings=dict(allow_extra_args=True),
)
@click.argument("benchmark_path", type=click.Path(exists=True))
@click.option(
    "-m",
    "--module",
    help="Module id to execute. If specified, runs a specific module instead of the full benchmark.",
    type=str,
    default=None,
)
@click.option(
    "-i",
    "--input-dir",
    help="Path to the folder with the appropriate input files (only used with --module).",
    type=click.Path(exists=True, writable=True),
    default=None,
)
@click.option(
    "-c",
    "--cores",
    help="Use at most N CPU cores in parallel. Default is 1.",
    type=int,
    default=1,
)
@click.option(
    "-u",
    "--update",
    help="Force re-run execution for all modules and stages.",
    is_flag=True,
    default=False,
)
@click.option("-d", "--dry", help="Dry run.", is_flag=True, default=False)
@click.option(
    "-y",
    "--yes",
    help="Automatically confirm all prompts.",
    is_flag=True,
    default=False,
)
@click.option(
    "-k",
    "--continue-on-error",
    help="Go on with independent jobs if a job fails (--keep-going in snakemake).",
    is_flag=True,
    default=False,
)
@click.option(
    "-r",
    "--use-remote-storage",
    help="Execute and store results remotely. Default False.",
    is_flag=True,
    default=False,
)
@click.option(
    "--keep-module-logs/--no-keep-module-logs",
    default=False,
    help="Keep module-specific log files after execution.",
)
@click.option(
    "--task-timeout",
    type=str,
    default=None,
    help="A `human friendly` timeout for each separate task execution (local only). Do note that total runtime is not additive. Example: 4h, 42m, 12s",
)
@click.option(
    "--executor",
    help="Specify a custom executor to use for workflow execution. Example: slurm",
    type=str,
    default=None,
)
@click.option(
    "--out-dir",
    help="Output folder name (local only). Default: `out`",
    default=None,
    type=str,
)
@click.pass_context
def run(
    ctx,
    benchmark_path,
    module,
    input_dir,
    cores,
    update,
    dry,
    yes,
    use_remote_storage,
    keep_module_logs,
    continue_on_error,
    task_timeout,
    executor,
    out_dir,
):
    """Run a benchmark or a specific module.

    Examples:

      ob run benchmark.yaml                    # Run full benchmark

      ob run benchmark.yaml --module my_module --input-dir ./data  # Run specific module
    """
    ctx.ensure_object(dict)
    extra_args = parse_extra_args(ctx.args)

    # Retrieve the global debug flag from the Click context
    debug = ctx.obj.get("DEBUG", False)

    # Parse timeout
    if task_timeout is None:
        timeout_s = None
    else:
        try:
            timeout_s = int(humanfriendly.parse_timespan(task_timeout))
        except humanfriendly.InvalidTimespan:
            logger.error(f"Invalid timeout value: {task_timeout}")
            sys.exit(1)

    # Decide whether to run whole benchmark or module
    if module:
        _run_module(
            benchmark_path=benchmark_path,
            module=module,
            input_dir=input_dir,
            out_dir=out_dir,
            out_dir_explicit=out_dir is not None,
            cores=cores,
            update=update,
            dry=dry,
            keep_module_logs=keep_module_logs,
            continue_on_error=continue_on_error,
            timeout_s=timeout_s,
            executor=executor,
            extra_args=extra_args,
        )
    else:
        # Validate out_dir usage
        if use_remote_storage and out_dir:
            raise click.UsageError("--out-dir can only be used with local storage")

        _run_benchmark(
            benchmark_path=benchmark_path,
            cores=cores,
            update=update,
            dry=dry,
            yes=yes,
            use_remote_storage=use_remote_storage,
            keep_module_logs=keep_module_logs,
            continue_on_error=continue_on_error,
            timeout_s=timeout_s,
            executor=executor,
            out_dir=out_dir,
            debug=debug,
            extra_args=extra_args,
        )


def _run_benchmark(
    benchmark_path,
    cores,
    update,
    dry,
    yes,
    use_remote_storage,
    keep_module_logs,
    continue_on_error,
    timeout_s,
    executor,
    out_dir,
    debug,
    extra_args,
):
    """Run a full benchmark."""
    try:
        out_dir = out_dir if out_dir else "out"
        b = BenchmarkExecution(Path(benchmark_path), Path(out_dir))
        logger.info("Benchmark YAML file integrity check passed.")
    except BenchmarkParseError as e:
        # Format parse errors with file location and context
        formatted_error = pretty_print_parse_error(e)
        log_error_and_quit(logger, f"Failed to load benchmark: {formatted_error}")
        return
    except Exception as e:
        log_error_and_quit(logger, f"Failed to load benchmark: {e}")
        return
    assert b is not None

    # it is as --continue_on_error originally
    if continue_on_error and not yes:
        msg = "run the full benchmark even if some jobs fail"
        abort_if_user_does_not_confirm(msg, logger)

    if update and not dry and not yes:
        msg = "re-run the entire workflow"
        abort_if_user_does_not_confirm(msg, logger)

    if use_remote_storage:
        storage_options = remote_storage_snakemake_args(b)
        # creates bucket if it doesn't exist
        _ = get_storage_from_benchmark(b)
    else:
        storage_options = {}

    workflow: WorkflowEngine = SnakemakeEngine()

    # Controlling resource allocation with Snakemake is tricky
    # -c only controls the number of parallelism for the Snakemake scheduler,
    # but there's no built-in mechanism in local execution to enforce well-behaved resource
    # allocation. Adding to that, bioinformatics tools will not always allow to control cores/threads

    logger.info("Running benchmark...")
    success = workflow.run_workflow(
        b,
        cores=cores,
        update=update,
        dryrun=dry,
        continue_on_error=continue_on_error,
        keep_module_logs=keep_module_logs,
        backend=b.get_benchmark_software_backend(),
        debug=debug,
        executor=executor,
        local_timeout=timeout_s,
        **storage_options,
        **extra_args,
    )

    log_result_and_quit(logger, success, type="Benchmark")


def _run_module(
    benchmark_path,
    module,
    input_dir,
    out_dir,
    out_dir_explicit,
    cores,
    update,
    dry,
    keep_module_logs,
    continue_on_error,
    timeout_s,
    executor,
    extra_args,
):
    """Run a specific module that is part of the benchmark."""
    logger.info("Running module on a local dataset.")

    # Use provided out_dir or default to "out"
    work_dir = Path(out_dir) if out_dir else Path("out")

    # Convert paths to absolute paths to avoid issues when work_dir changes
    work_dir = work_dir.resolve()
    benchmark_path_abs = Path(benchmark_path).resolve()

    try:
        b = BenchmarkExecution(benchmark_path_abs, work_dir)
    except BenchmarkParseError as e:
        formatted_error = pretty_print_parse_error(e)
        log_error_and_quit(logger, f"Failed to load benchmark: {formatted_error}")
        return
    except Exception as e:
        log_error_and_quit(logger, f"Failed to load benchmark: {e}")
        return
    assert b is not None

    benchmark_nodes = b.get_nodes_by_module_id(module_id=module)

    is_entrypoint_module = all([node.is_entrypoint() for node in benchmark_nodes])

    if len(benchmark_nodes) <= 0:
        log_error_and_quit(
            logger,
            f"Error: Could not find module with id `{module}` in benchmark definition",
        )
        return

    assert len(benchmark_nodes) > 0
    logger.info(f"Found {len(benchmark_nodes)} workflow nodes for module {module}.")

    if input_dir is None:
        # Default to 'out' folder in current directory if it exists
        default_input_dir = os.path.join(os.getcwd(), "out")
        if os.path.exists(default_input_dir) and os.path.isdir(default_input_dir):
            input_dir = default_input_dir
            logger.info(f"Using default input directory: {input_dir}")
        elif not is_entrypoint_module:
            log_error_and_quit(
                logger,
                "Error: --input-dir is required for non-entrypoint modules when 'out' folder doesn't exist in current directory.",
            )
            return
        else:
            input_dir = os.getcwd()

    logger.info("Running module benchmark...")

    # Check if input path exists and is a directory
    if not (os.path.exists(input_dir) and os.path.isdir(input_dir)):
        log_error_and_quit(
            logger,
            f"Error: Input directory does not exist or is not a valid directory: `{input_dir}`",
        )
        return
    assert os.path.exists(input_dir) and os.path.isdir(input_dir)

    benchmark_datasets = b.get_benchmark_datasets()

    # Initialize actual_input_dir to track where files are actually located
    actual_input_dir = input_dir

    # Check available files in input to figure out what dataset are we processing
    if module in benchmark_datasets:
        # we're given the initial dataset module to process, then we know
        dataset = module
    else:
        # we try to figure the dataset based on the files present in the input directory
        files = os.listdir(actual_input_dir)
        base_names = [file.split(".")[0] for file in files]
        dataset = next((d for d in benchmark_datasets if d in base_names), None)

        # If not found directly, look in data/ subdirectory structure: data/{dataset}/{params}/
        if dataset is None:
            data_dir = os.path.join(input_dir, "data")
            if os.path.exists(data_dir) and os.path.isdir(data_dir):
                subdirs = os.listdir(data_dir)
                dataset = next((d for d in benchmark_datasets if d in subdirs), None)
                if dataset is not None:
                    # Update actual_input_dir to point to the actual data location
                    # Find the first parameter directory under data/{dataset}/
                    dataset_dir = os.path.join(data_dir, dataset)
                    if os.path.exists(dataset_dir) and os.path.isdir(dataset_dir):
                        param_dirs = os.listdir(dataset_dir)
                        if len(param_dirs) > 0:
                            actual_input_dir = os.path.join(dataset_dir, param_dirs[0])

    if dataset is None:
        log_error_and_quit(
            logger,
            f"Error: Could not infer the appropriate dataset to run the node workflow on based on the files available in `{input_dir}`. None of the available datasets {benchmark_datasets} match the base names of the files.",
        )
        return

    assert dataset is not None

    # Check if input directory contains all necessary input files
    required_inputs = list(map(lambda node: node.get_inputs(), benchmark_nodes))
    required_inputs = list(chain.from_iterable(required_inputs))
    required_input_files = list(
        set([os.path.basename(path) for path in required_inputs])
    )
    required_input_files = [
        file.format(dataset=dataset) for file in required_input_files
    ]

    input_files = os.listdir(actual_input_dir)
    missing_files = [file for file in required_input_files if file not in input_files]

    if len(missing_files) > 0:
        log_error_and_quit(
            logger,
            f"Error: The following required input files are missing from the input directory: {missing_files}.",
        )
        return
    assert len(missing_files) == 0

    # If out_dir was explicitly provided and differs from input location,
    # create output structure and symlink input files
    # This allows Snakemake to find inputs while writing outputs to the correct location
    workflow_input_dir = actual_input_dir
    if (
        out_dir_explicit
        and Path(actual_input_dir).resolve() != (work_dir / "data" / dataset).resolve()
    ):
        # Determine the relative path structure from actual_input_dir
        # Expected: actual_input_dir is something like input_dir/data/{dataset}/{params}
        try:
            # Extract the path components
            input_path = Path(actual_input_dir).resolve()
            # Find the data/ directory level
            parts = input_path.parts
            data_idx = next((i for i, p in enumerate(parts) if p == "data"), None)

            if data_idx is not None and data_idx + 2 < len(parts):
                # Reconstruct: data/{dataset}/{params}
                relative_structure = Path(*parts[data_idx:])
                output_data_dir = work_dir / relative_structure

                # Create the output directory structure
                os.makedirs(output_data_dir, exist_ok=True)

                # Copy only the required input FILES (not directories) to the output directory
                # This ensures we don't copy output directories from previous module runs
                # Using copy instead of symlink avoids issues with container bind mounts
                for file in required_input_files:
                    src = Path(actual_input_dir) / file
                    dst = output_data_dir / file
                    # Only copy if it's a file (not a directory) and doesn't already exist
                    if src.is_file() and not dst.exists():
                        shutil.copy2(src, dst)
                        logger.info(f"Copied {src} -> {dst}")

                # Recreate human-readable parameter symlinks in the dataset directory
                # These are symlinks like dataset_generator-fcps_dataset_name-atom -> .46bb23e0
                dataset_dir = output_data_dir.parent  # This is work_dir/data/{dataset}
                source_dataset_dir = Path(
                    actual_input_dir
                ).parent  # Original data/{dataset} dir
                if source_dataset_dir.exists():
                    for item in source_dataset_dir.iterdir():
                        if item.is_symlink():
                            # Recreate the symlink in the new location
                            target = item.readlink()
                            new_symlink = dataset_dir / item.name
                            if not new_symlink.exists():
                                os.symlink(target, new_symlink)
                                logger.info(
                                    f"Recreated symlink {new_symlink} -> {target}"
                                )

                # Update workflow_input_dir to point to the copied location
                workflow_input_dir = str(output_data_dir)
        except Exception as e:
            logger.warning(
                f"Could not create copied input structure: {e}. Using original input directory."
            )

    logger.info(f"Running workflow with input_dir: {workflow_input_dir}")
    logger.info(f"Dataset: {dataset}")
    logger.info(f"Work dir: {work_dir if out_dir_explicit else 'default (cwd)'}")

    workflow: WorkflowEngine = SnakemakeEngine()

    # Prepare workflow arguments (without node, which is set in the loop)
    workflow_args = {
        "input_dir": Path(workflow_input_dir),
        "dataset": dataset,
        "cores": 1,
        "update": update,
        "dryrun": dry,
        "continue_on_error": continue_on_error,
        "keep_module_logs": keep_module_logs,
        "backend": b.get_benchmark_software_backend(),
        "executor": executor,
        "local_timeout": timeout_s,
        "benchmark_file_path": b.get_definition_file(),
    }

    # Only pass work_dir if out_dir was explicitly provided
    if out_dir_explicit:
        workflow_args["work_dir"] = work_dir

    workflow_args.update(extra_args)

    for benchmark_node in benchmark_nodes:
        # When running a single module, it doesn't have sense to make parallelism level (cores) configurable
        workflow_args["node"] = benchmark_node
        success = workflow.run_node_workflow(**workflow_args)

        log_result_and_quit(logger, success, type="Module")


def log_error_and_quit(logger, error):
    logger.error(error)
    sys.exit(1)


def log_result_and_quit(logger, success: bool, type: str):
    if success:
        logger.info(f"{type} run has finished successfully.")
        sys.exit(0)
    logger.error(f"{type} run has failed.")
    sys.exit(1)


def abort_if_user_does_not_confirm(msg: str, logger):
    _msg = f"Are you sure you want to {msg}?"
    if not click.confirm(_msg, abort=True):
        logger.debug("aborting")
        raise click.Abort()
