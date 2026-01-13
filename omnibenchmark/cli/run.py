"""cli commands related to benchmark/module execution and start"""

import os
import shutil
import sys
from itertools import chain
from pathlib import Path

import click
import humanfriendly
from pydantic import ValidationError as PydanticValidationError

from omnibenchmark.backend.collector_resolution import resolve_metric_collectors
from omnibenchmark.backend.resolver import ModuleResolver
from omnibenchmark.backend.snakemake_gen import SnakemakeGenerator, save_metadata
from omnibenchmark.benchmark import BenchmarkExecution
from omnibenchmark.benchmark.params import Params
from omnibenchmark.benchmark.symlinks import SymlinkManager
from omnibenchmark.cli.error_formatting import pretty_print_parse_error
from omnibenchmark.cli.utils.args import parse_extra_args
from omnibenchmark.cli.utils.logging import logger
from omnibenchmark.config import get_git_cache_dir
from omnibenchmark.git import get_or_update_cached_repo
from omnibenchmark.model.resolved import ResolvedNode
from omnibenchmark.model.validation import BenchmarkParseError
from omnibenchmark.remote.storage import (
    get_storage_from_benchmark,
    remote_storage_snakemake_args,
)
from omnibenchmark.workflow.snakemake import SnakemakeEngine
from omnibenchmark.workflow.workflow import WorkflowEngine


def format_pydantic_errors(e: PydanticValidationError) -> str:
    """Format Pydantic validation errors to show which fields are missing or invalid."""
    error_lines = ["Validation failed:"]
    for error in e.errors():
        field = " -> ".join(str(loc) for loc in error["loc"])
        msg = error["msg"]
        error_type = error["type"]

        # Make missing field errors more explicit
        if error_type == "missing":
            error_lines.append(f"  - Missing required field: '{field}'")
        else:
            error_lines.append(f"  - Field '{field}': {msg}")

    return "\n".join(error_lines)


@click.command(
    name="run",
    context_settings=dict(allow_extra_args=True),
)
@click.argument("benchmark", type=click.Path(exists=True))
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
@click.option(
    "--only-snakefile",
    is_flag=True,
    help="Generate executable Snakefile without running it (requires --dry). Produces runnable workflow for manual execution.",
)
@click.pass_context
def run(
    ctx,
    benchmark,
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
    only_snakefile,
):
    """Run a benchmark or a specific module.

    BENCHMARK: Path to benchmark YAML file.

    Examples:

      ob run benchmark.yaml                    # Run full benchmark

      ob run benchmark.yaml --module my_module --input-dir ./data  # Run specific module
    """
    ctx.ensure_object(dict)
    extra_args = parse_extra_args(ctx.args)

    # Retrieve the global debug flag from the Click context
    debug = ctx.obj.get("DEBUG", False)

    # Validate flag combinations
    if only_snakefile and not dry:
        raise click.UsageError(
            "--only-snakefile requires --dry flag (used to generate Snakefile without execution)"
        )

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
            benchmark_path=benchmark,
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
            benchmark_path=benchmark,
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
            only_snakefile=only_snakefile,
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
    only_snakefile,
):
    """Run a full benchmark."""
    try:
        out_dir = out_dir if out_dir else "out"
        benchmark_path_abs = Path(benchmark_path).resolve()
        b = BenchmarkExecution(benchmark_path_abs, Path(out_dir))
        logger.info("Benchmark YAML file integrity check passed.")
    except PydanticValidationError as e:
        # Handle Pydantic validation errors with detailed field information
        log_error_and_quit(
            logger, f"Failed to load benchmark:\n{format_pydantic_errors(e)}"
        )
        return
    except BenchmarkParseError as e:
        # Format parse errors with file location and context
        formatted_error = pretty_print_parse_error(e)
        log_error_and_quit(logger, f"Failed to load benchmark: {formatted_error}")
        return
    except Exception as e:
        log_error_and_quit(logger, f"Failed to load benchmark: {e}")
        return
    assert b is not None

    # NEW EXECUTION PATH: Generate explicit Snakefile during dry-run
    # This does:
    # 1. Populate git cache (fetch all repos)
    # 2. Clone modules to work directory (out/.repos/)
    # 3. Dereference entrypoints
    # 4. Generate explicit Snakefile
    # 5. Save metadata
    if dry:
        # Use clean UI when not in debug mode (debug flag is the opposite of clean UI)
        use_clean_ui = not debug

        # Start timing for the entire generation process
        import time

        start_time = time.time()

        _populate_git_cache(b, quiet=use_clean_ui, cores=cores)
        _generate_explicit_snakefile(
            benchmark=b,
            benchmark_yaml_path=benchmark_path_abs,
            out_dir=Path(out_dir),
            debug_mode=not only_snakefile,  # If --only-snakefile, generate executable; otherwise debug
            cores=cores,  # Pass cores for parallel module resolution
            quiet=use_clean_ui,  # Use clean progress UI when not debugging
            start_time=start_time,  # Pass start time for final report
        )

        # If --only-snakefile, exit without running
        if only_snakefile:
            if not use_clean_ui:
                logger.info("\nExecutable Snakefile generation complete. Exiting.")
                logger.info(
                    f"To execute: cd {out_dir} && snakemake --use-conda --cores <N>"
                )
            sys.exit(0)

        # Otherwise continue to run the generated debug Snakefile
        logger.info("\nDebug Snakefile generated. Now executing...")

    # Validation for interactive flags
    if continue_on_error and not yes and not dry:
        msg = "run the full benchmark even if some jobs fail"
        abort_if_user_does_not_confirm(msg, logger)

    if update and not yes and not dry:
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
    except PydanticValidationError as e:
        # Handle Pydantic validation errors with detailed field information
        log_error_and_quit(
            logger, f"Failed to load benchmark:\n{format_pydantic_errors(e)}"
        )
        return
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


def _populate_git_cache(
    benchmark: BenchmarkExecution, quiet: bool = False, cores: int = 4
):
    """Populate git cache with repositories from benchmark.

    This is part of the NEW execution path using the two-tier caching system.
    During dry-run, we populate the cache so subsequent real runs can use it.

    Args:
        benchmark: BenchmarkExecution instance
        quiet: If True, use progress bar instead of logging
        cores: Number of parallel workers for fetching
    """
    from omnibenchmark.cli.progress import ProgressDisplay

    cache_dir = get_git_cache_dir()

    # Extract unique repositories with their commits (for optimization)
    # Format: {repo_url: commit}
    repos = {}
    for stage in benchmark.model.stages:
        for module in stage.modules:
            if hasattr(module, "repository") and module.repository:
                repo_url = module.repository.url
                commit = module.repository.commit
                if repo_url:
                    repos[repo_url] = commit

    # Also include metric collector repositories
    if (
        hasattr(benchmark.model, "metric_collectors")
        and benchmark.model.metric_collectors
    ):
        for collector in benchmark.model.metric_collectors:
            if hasattr(collector, "repository") and collector.repository:
                repo_url = collector.repository.url
                commit = collector.repository.commit
                if repo_url:
                    repos[repo_url] = commit

    if not repos:
        if not quiet:
            logger.info("No repositories found in benchmark definition")
        return

    if not quiet:
        logger.info(f"Populating cache with {len(repos)} repositories...")

    progress = ProgressDisplay() if quiet else None
    if quiet:
        progress.start_task("Fetching repositories to cache", total=len(repos))

    # Populate cache in parallel
    success_count = 0
    failed = []

    from concurrent.futures import ThreadPoolExecutor, as_completed
    import threading

    lock = threading.Lock()

    def fetch_repo_task(repo_url, commit):
        """Fetch a single repository (runs in parallel)."""
        # Check if we need to fetch at all
        from omnibenchmark.git.cache import parse_repo_url

        repo_path = parse_repo_url(repo_url)
        repo_cache_dir = cache_dir / repo_path

        # Skip fetch if:
        # 1. Repo exists in cache
        # 2. Commit looks like a commit hash (40 hex chars)
        # 3. We have that commit locally
        skip_fetch = False
        if (
            repo_cache_dir.exists()
            and commit
            and len(commit) == 40
            and all(c in "0123456789abcdef" for c in commit.lower())
        ):
            try:
                from dulwich import porcelain

                repo = porcelain.open_repo(str(repo_cache_dir))
                # Check if we have this commit
                repo[commit.encode("ascii")]
                skip_fetch = True
            except (KeyError, Exception):
                pass

        if skip_fetch:
            return (repo_url, True, None)

        if not quiet:
            logger.info(f"Caching {repo_url}")

        # Suppress dulwich logging in quiet mode (thread-safe approach)
        if quiet:
            import logging

            dulwich_logger = logging.getLogger("dulwich")
            old_level = dulwich_logger.level
            dulwich_logger.setLevel(logging.ERROR)

        try:
            get_or_update_cached_repo(repo_url, cache_dir)
            return (repo_url, True, None)
        except Exception as e:
            return (repo_url, False, str(e))
        finally:
            if quiet:
                dulwich_logger.setLevel(old_level)

    with ThreadPoolExecutor(max_workers=cores) as executor:
        futures = {
            executor.submit(fetch_repo_task, repo_url, commit): repo_url
            for repo_url, commit in repos.items()
        }

        for future in as_completed(futures):
            repo_url = futures[future]

            try:
                repo_url_result, success, error = future.result()

                if success:
                    with lock:
                        success_count += 1
                else:
                    logger.error(f"Failed to cache {repo_url}: {error}")
                    with lock:
                        failed.append((repo_url, error))
            except Exception as e:
                # Catch any unhandled exceptions from the task
                logger.error(f"Exception while caching {repo_url}: {e}")
                with lock:
                    failed.append((repo_url, str(e)))

            if quiet:
                progress.update(advance=1)

    if quiet:
        progress.finish()
        progress.success(f"Cached {success_count}/{len(repos)} repositories")
    else:
        logger.info(f"Successfully cached {success_count}/{len(repos)} repositories")

    if failed:
        logger.warning(f"Failed to cache {len(failed)} repositories:")
        for repo_url, error in failed:
            logger.warning(f"  {repo_url}: {error}")


def _prepare_parameter_directories(resolved_nodes, out_dir, quiet: bool = False):
    """
    Pre-create parameter directories and human-readable symlinks.

    For each parameterized node, this:
    1. Creates the .{hash} directory for content-addressable storage
    2. Creates a human-readable symlink pointing to the .{hash} directory
    3. Writes parameters.json inside the .{hash} directory

    Args:
        resolved_nodes: List of ResolvedNode instances
        out_dir: Output directory base path
        quiet: If True, suppress logging
    """
    if not quiet:
        logger.info("\nPre-creating parameter directories and symlinks:")

    created_count = 0
    for node in resolved_nodes:
        if node.parameters:
            # Determine base directory for this node's parameter directory
            # The node's param_id (like .515ef2c0) tells us which hash directory is ours
            # We need to find where in the output path our param_id appears
            if node.outputs:
                output_path = Path(node.outputs[0])
                parts = output_path.parts

                # Find our param_id in the path
                param_id = node.param_id  # e.g., ".515ef2c0"
                try:
                    param_idx = parts.index(param_id)
                    # Base directory is everything before our param_id
                    base_dir = out_dir / Path(*parts[:param_idx])
                except ValueError:
                    # param_id not found in path - this shouldn't happen
                    logger.warning(
                        f"  Could not find param_id {param_id} in output path {output_path}"
                    )
                    continue

                # Create symlink manager for this base directory
                manager = SymlinkManager(base_dir=base_dir)

                # Store parameters (creates .hash dir and human-readable symlink)
                try:
                    symlink_info = manager.store(node.parameters)
                    logger.debug(
                        f"  {node.stage_id}/{node.module_id}: {symlink_info['human']} -> {symlink_info['folder']}"
                    )
                    created_count += 1
                except Exception as e:
                    logger.warning(f"  Failed to create symlink for {node.id}: {e}")

    if not quiet:
        logger.info(
            f"Pre-created directories and symlinks for {created_count} parameterized nodes"
        )


def _generate_explicit_snakefile(
    benchmark: BenchmarkExecution,
    benchmark_yaml_path: Path,
    out_dir: Path,
    nesting_strategy: str = "nested",
    debug_mode: bool = True,
    cores: int = 4,
    quiet: bool = False,
    start_time: float = None,
):
    """
    Generate explicit Snakefile during dry-run.

    This resolves all modules, dereferences entrypoints, and generates
    a human-readable Snakefile with shell: directives.

    Args:
        benchmark: BenchmarkExecution instance
        benchmark_yaml_path: Path to original benchmark YAML
        out_dir: Output directory
        nesting_strategy: Path nesting strategy - "nested" (default) or "flat"
        cores: Number of parallel workers for module resolution
        quiet: If True, use clean progress UI; if False, show verbose INFO logs
        start_time: Start time for timing (if None, starts now)
    """
    from omnibenchmark.cli.progress import ProgressDisplay
    import time

    progress = ProgressDisplay()
    if start_time is None:
        start_time = time.time()

    if not quiet:
        logger.info("\nGenerating explicit Snakefile...")

    # Create resolver with work dir relative to output
    work_dir = out_dir / ".modules"

    # Get benchmark directory (where benchmark.yaml lives)
    benchmark_dir = benchmark_yaml_path.parent

    # Initialize resolver with software environment information
    resolver = ModuleResolver(
        work_base_dir=work_dir,
        output_dir=out_dir,
        software_backend=benchmark.model.get_software_backend(),
        software_environments=benchmark.model.get_software_environments(),
        benchmark_dir=benchmark_dir,
    )

    # Collect all unique modules to resolve
    unique_modules = {}  # module_id -> (module, software_environment_id)
    for stage in benchmark.model.stages:
        for module in stage.modules:
            if module.id not in unique_modules:
                unique_modules[module.id] = (module, module.software_environment)

    # Resolve modules in parallel with progress bar
    if not quiet:
        logger.info(f"\nResolving {len(unique_modules)} modules...")

    from concurrent.futures import ThreadPoolExecutor, as_completed
    import threading

    resolved_modules_cache = {}
    resolution_lock = threading.Lock()
    captured_warnings = []

    def resolve_module_task(module_id, module, software_env_id):
        """Resolve a single module (runs in parallel)."""
        import warnings
        import logging

        # Only suppress logging if in quiet mode (don't touch stdout/stderr as it breaks Rich)
        if quiet:
            # Disable all logging by removing all handlers temporarily
            root_logger = logging.getLogger()
            old_handlers = root_logger.handlers[:]
            old_level = root_logger.level
            root_logger.handlers = []
            root_logger.setLevel(logging.CRITICAL + 1)

            # Disable all loggers in the omnibenchmark namespace
            omnibenchmark_logger = logging.getLogger("omnibenchmark")
            old_omni_level = omnibenchmark_logger.level
            old_omni_propagate = omnibenchmark_logger.propagate
            old_omni_handlers = omnibenchmark_logger.handlers[:]
            omnibenchmark_logger.handlers = []
            omnibenchmark_logger.setLevel(logging.CRITICAL + 1)
            omnibenchmark_logger.propagate = False

        # Capture warnings instead of suppressing
        module_warnings = []
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            if quiet:
                # Suppress warning display (they're still captured in w)
                warnings.filterwarnings("ignore")

            try:
                resolved = resolver.resolve(
                    module=module,
                    module_id=module_id,
                    software_environment_id=software_env_id,
                    dirty=False,
                )
                # Capture any warnings
                for warning in w:
                    module_warnings.append((module_id, str(warning.message)))

                return (module_id, resolved, None, module_warnings)
            except Exception as e:
                return (module_id, None, str(e), module_warnings)
            finally:
                # Restore logging if it was suppressed
                if quiet:
                    root_logger.handlers = old_handlers
                    root_logger.setLevel(old_level)
                    omnibenchmark_logger.handlers = old_omni_handlers
                    omnibenchmark_logger.setLevel(old_omni_level)
                    omnibenchmark_logger.propagate = old_omni_propagate

    if quiet:
        progress.start_task("Resolving modules", total=len(unique_modules))

    with ThreadPoolExecutor(max_workers=cores) as executor:
        # Submit all resolution tasks
        futures = {
            executor.submit(resolve_module_task, mod_id, mod, env_id): mod_id
            for mod_id, (mod, env_id) in unique_modules.items()
        }

        # Process completed tasks
        for future in as_completed(futures):
            module_id, resolved, error, module_warnings = future.result()

            # Collect warnings
            with resolution_lock:
                captured_warnings.extend(module_warnings)

            if error:
                if quiet:
                    progress.error(f"{module_id}: {error}")
                else:
                    logger.error(f"Failed to resolve {module_id}: {error}")
            else:
                with resolution_lock:
                    resolved_modules_cache[module_id] = resolved

            if quiet:
                progress.update(advance=1)

    if quiet:
        progress.finish()

        if len(resolved_modules_cache) < len(unique_modules):
            progress.error(
                f"Failed to resolve {len(unique_modules) - len(resolved_modules_cache)} modules"
            )
        else:
            progress.success(f"Resolved {len(resolved_modules_cache)} modules")

        # Display captured warnings at the end
        if captured_warnings:
            progress.console.print()
            progress.console.print(
                f"[yellow]⚠ {len(captured_warnings)} warning(s) during resolution:[/yellow]"
            )
            for module_id, warning_msg in captured_warnings:
                progress.console.print(
                    f"  [yellow]•[/yellow] [dim]{module_id}:[/dim] {warning_msg}"
                )
    else:
        logger.info(
            f"Successfully resolved {len(resolved_modules_cache)}/{len(unique_modules)} modules"
        )

    # Build nodes from resolved modules
    if not quiet:
        logger.info("\nBuilding execution graph...")

    # Collect resolved nodes
    resolved_nodes = []

    # Create a mapping of output_id -> list of (node_id, output_path) for dependency resolution
    output_to_nodes = {}

    # Track nodes from previous stage for cartesian product
    previous_stage_nodes = []

    for stage in benchmark.model.stages:
        # Nodes created in this stage (for next stage's cartesian product)
        current_stage_nodes = []

        for module in stage.modules:
            module_id = module.id

            try:
                # Get resolved module from cache
                if module_id not in resolved_modules_cache:
                    logger.warning(f"      Module {module_id} not in cache, skipping")
                    continue

                resolved_module = resolved_modules_cache[module_id]

                # Expand parameters
                if module.parameters:
                    params_list = []
                    for param in module.parameters:
                        params_list.extend(Params.expand_from_parameter(param))
                else:
                    params_list = [None]

                # Get input combinations from previous stage
                # If this stage has inputs, create cartesian product with previous stage nodes
                # If no inputs (initial stage), just iterate over our own parameters
                if stage.inputs and previous_stage_nodes:
                    # Cartesian product: each previous node × each param combination
                    from itertools import product

                    input_node_combinations = previous_stage_nodes
                    node_combinations = list(
                        product(input_node_combinations, params_list)
                    )
                else:
                    # Initial stage: just our parameter combinations
                    node_combinations = [(None, params) for params in params_list]

                # Create a node for each combination
                for input_node, params in node_combinations:
                    # Check if this combination should be excluded
                    if input_node and module.exclude:
                        # Check if input node's module is in the exclude list
                        if input_node.module_id in module.exclude:
                            logger.debug(
                                f"      Excluding combination: {input_node.module_id} → {module_id}"
                            )
                            continue

                    param_id = f".{params.hash_short()}" if params else ".default"

                    # Build node ID
                    if input_node:
                        # Include input node ID in our ID for uniqueness
                        node_id = f"{input_node.id}-{stage.id}-{module_id}{param_id}"
                    else:
                        node_id = f"{stage.id}-{module_id}{param_id}"

                    # Resolve inputs from the input_node
                    inputs = {}
                    input_name_mapping = {}  # Map sanitized -> original names
                    base_path = None

                    if input_node:
                        # Build a mapping of output_id -> output_path for this node and its ancestors
                        # This allows us to resolve inputs from any stage in the DAG
                        output_id_to_path = {}

                        # Helper to extract output_id from output path
                        def get_output_ids_for_node(node):
                            """Get output_id -> path mapping for a node by matching against stage outputs."""
                            result = {}
                            # Find the stage this node belongs to
                            for s in benchmark.model.stages:
                                if s.id == node.stage_id:
                                    # Match node's output paths to stage's output specs
                                    for output_spec in s.outputs:
                                        for output_path in node.outputs:
                                            # Check if this output_path matches the output_spec pattern
                                            # We can check if the spec's path template appears in the output_path
                                            # For now, use a simple index-based mapping
                                            if s.outputs.index(output_spec) < len(
                                                node.outputs
                                            ):
                                                result[output_spec.id] = node.outputs[
                                                    s.outputs.index(output_spec)
                                                ]
                            return result

                        # Collect outputs from input_node
                        output_id_to_path.update(get_output_ids_for_node(input_node))

                        # Traverse back through the lineage to collect ancestor outputs
                        # Extract the lineage by parsing the node_id which contains parent IDs
                        # For a node like "data_clustbench_46bb23e0_clustering_fastcluster_...",
                        # the lineage is encoded in the ID prefix
                        def get_ancestor_nodes(node):
                            """Get all ancestor nodes by tracing back through the lineage."""
                            ancestors = []
                            # Parse node ID to find parent nodes
                            # Node IDs are like: parent_id + "-" + stage_id + "-" + module_id + param_id
                            # We need to find all nodes whose ID is a prefix of this node's ID
                            for prev_node in resolved_nodes:
                                if node.id.startswith(prev_node.id + "-"):
                                    ancestors.append(prev_node)
                                    # Recursively get ancestors of this ancestor
                                    ancestors.extend(get_ancestor_nodes(prev_node))
                            return ancestors

                        # Collect outputs from all ancestors
                        for ancestor in get_ancestor_nodes(input_node):
                            output_id_to_path.update(get_output_ids_for_node(ancestor))

                        # Now map stage inputs to actual paths using output_id
                        if stage.inputs:
                            for input_collection in stage.inputs:
                                for input_id in input_collection.entries:
                                    # input_id is something like "data.true_labels" or "clustering.predicted_ks_range"
                                    # Try to find a matching output_id
                                    if input_id in output_id_to_path:
                                        sanitized_id = input_id.replace(".", "_")
                                        inputs[sanitized_id] = output_id_to_path[
                                            input_id
                                        ]
                                        input_name_mapping[sanitized_id] = input_id
                                    else:
                                        logger.warning(
                                            f"      Could not resolve input {input_id} for node {node_id}"
                                        )

                        # Use the first input's directory as base for nesting
                        if inputs:
                            from pathlib import Path as PathLib

                            first_input = next(iter(inputs.values()))
                            base_path = str(PathLib(first_input).parent)

                    # Determine dataset value BEFORE building outputs
                    # For entrypoint nodes (no inputs), use module_id as dataset
                    # For subsequent nodes, propagate from input_node if available
                    if not inputs:
                        # Entrypoint node - dataset is the module_id
                        dataset_value = module_id
                    elif (
                        input_node
                        and hasattr(input_node, "dataset")
                        and input_node.dataset
                    ):
                        # Propagate dataset from input node
                        dataset_value = input_node.dataset
                    else:
                        # Fallback - no dataset value available
                        dataset_value = None

                    # Build output paths based on nesting strategy
                    outputs = []
                    for output_spec in stage.outputs:
                        # Get the output path template
                        output_path_template = output_spec.path

                        # Replace {dataset} placeholder with actual dataset value
                        # This is NOT a Snakemake wildcard - it's an omnibenchmark path template variable
                        if dataset_value and "{dataset}" in output_path_template:
                            output_path_template = output_path_template.replace(
                                "{dataset}", dataset_value
                            )

                        if nesting_strategy == "nested":
                            # Nested strategy: dependent stages nest under their input directories
                            if base_path:
                                # Has inputs: nest under input directory
                                output_path = f"{base_path}/{stage.id}/{module_id}/{param_id}/{output_path_template}"
                            else:
                                # Initial stage: stage/module/param/filename
                                output_path = f"{stage.id}/{module_id}/{param_id}/{output_path_template}"
                        elif nesting_strategy == "flat":
                            # Flat strategy: all stages at same level (stage/module/param/filename)
                            output_path = f"{stage.id}/{module_id}/{param_id}/{output_path_template}"
                        else:
                            raise ValueError(
                                f"Unknown nesting strategy: {nesting_strategy}"
                            )

                        outputs.append(output_path)
                        # Map output_id to this node for dependency resolution
                        if output_spec.id not in output_to_nodes:
                            output_to_nodes[output_spec.id] = []
                        output_to_nodes[output_spec.id].append((node_id, output_path))

                    # Determine resource requirements (module-level overrides stage-level)
                    # Hierarchy: module.resources > stage.resources > default (2 cores)
                    node_resources = None
                    if hasattr(module, "resources") and module.resources:
                        # Module-level resources take precedence
                        node_resources = module.resources
                    elif hasattr(stage, "resources") and stage.resources:
                        # Fall back to stage-level resources
                        node_resources = stage.resources
                    # If neither is set, node_resources remains None (will use default in Snakefile gen)

                    # Create ResolvedNode
                    node = ResolvedNode(
                        id=node_id,
                        stage_id=stage.id,
                        module_id=module_id,
                        param_id=param_id,
                        module=resolved_module,
                        parameters=params,
                        inputs=inputs,
                        outputs=outputs,
                        dataset=dataset_value,
                        input_name_mapping=input_name_mapping,
                        benchmark_name=benchmark.model.get_name(),
                        benchmark_version=benchmark.model.get_version(),
                        benchmark_author=benchmark.model.get_author(),
                        resources=node_resources,
                    )

                    resolved_nodes.append(node)
                    current_stage_nodes.append(node)

                if not quiet:
                    logger.info(
                        f"      Created {len([n for n in current_stage_nodes if n.module_id == module_id])} nodes for {module_id}"
                    )

            except Exception as e:
                logger.error(f"      Failed to resolve {module_id}: {e}")
                import traceback

                if logger.level <= 10:  # DEBUG level
                    traceback.print_exc()

        # Update previous_stage_nodes for next iteration
        previous_stage_nodes = current_stage_nodes

    if not quiet:
        logger.info(f"Created {len(resolved_nodes)} nodes")

    # Note: Parameter directories and JSON files are now created by the Snakefile itself
    # during execution, making the workflow self-contained and portable

    # ========== Resolve metric collectors as regular nodes ==========
    # Metric collectors are converted to ResolvedNode instances and added to the DAG
    collector_nodes = resolve_metric_collectors(
        metric_collectors=benchmark.model.get_metric_collectors(),
        resolved_nodes=resolved_nodes,
        benchmark=benchmark.model,
        resolver=resolver,
        quiet=quiet,
    )

    # Add collector nodes to the main resolved_nodes list
    # They will be treated as regular nodes in the Snakefile generation
    resolved_nodes.extend(collector_nodes)

    # Generate Snakefile (no status message needed)

    snakefile_path = out_dir / "Snakefile"

    generator = SnakemakeGenerator(
        benchmark_name=benchmark.model.get_name(),
        benchmark_version=benchmark.model.get_version(),
        benchmark_author=benchmark.model.get_author(),
    )

    generator.generate_snakefile(
        nodes=resolved_nodes,
        collectors=[],  # Empty - collectors are now in resolved_nodes
        output_path=snakefile_path,
        debug_mode=debug_mode,
    )

    # Save metadata
    save_metadata(
        benchmark_yaml_path=benchmark_yaml_path,
        output_dir=out_dir,
        nodes=resolved_nodes,
        collectors=[],  # Empty - collectors are now in resolved_nodes
    )

    if quiet:
        # Simple one-line summary with timing
        elapsed = time.time() - start_time
        progress.success(
            f"Generated {len(resolved_nodes)} rules in {elapsed:.1f}s in {snakefile_path}"
        )
    else:
        logger.info(f"\nSnakefile generation complete: {snakefile_path}")
        logger.info(f"  Benchmark: {benchmark.model.get_name()}")
        logger.info(f"  Modules: {len(resolved_modules_cache)}")
        logger.info(f"  Nodes: {len(resolved_nodes)}")
