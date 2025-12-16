"""cli commands related to benchmark/module execution and start"""

import sys
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
    """Run a benchmark as specified in the yaml."""
    ctx.ensure_object(dict)
    extra_args = parse_extra_args(ctx.args)

    # Retrieve the global debug flag from the Click context
    debug = ctx.obj.get("DEBUG", False)

    if task_timeout is None:
        timeout_s = None
    else:
        try:
            timeout_s = int(humanfriendly.parse_timespan(task_timeout))
        except humanfriendly.InvalidTimespan:
            logger.error(f"Invalid timeout value: {task_timeout}")
            sys.exit(1)

    # Validate out_dir usage
    if use_remote_storage and out_dir:
        raise click.UsageError("--out-dir can only be used with local storage")

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
