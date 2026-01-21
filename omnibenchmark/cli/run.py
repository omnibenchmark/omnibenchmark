"""cli commands related to benchmark/module execution and start"""

import os
import subprocess
import sys
from pathlib import Path

import click
from pydantic import ValidationError as PydanticValidationError

from omnibenchmark.backend.collector_resolution import resolve_metric_collectors
from omnibenchmark.backend.resolver import ModuleResolver
from omnibenchmark.backend.snakemake_gen import SnakemakeGenerator, save_metadata
from omnibenchmark.benchmark import BenchmarkExecution
from omnibenchmark.benchmark.params import Params
from omnibenchmark.cli.error_formatting import pretty_print_parse_error
from omnibenchmark.cli.utils.logging import logger
from omnibenchmark.config import get_git_cache_dir
from omnibenchmark.git import get_or_update_cached_repo
from omnibenchmark.model import SoftwareBackendEnum
from omnibenchmark.model.resolved import ResolvedNode
from omnibenchmark.model.validation import BenchmarkParseError


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
)
@click.argument("benchmark", type=click.Path(exists=True))
@click.option(
    "-c",
    "--cores",
    help="Use at most N CPU cores in parallel. Default is 1.",
    type=int,
    default=1,
)
@click.option(
    "-d",
    "--dry",
    help="Dry run (only generate Snakefile, don't execute).",
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
    "--out-dir",
    help="Output folder name. Default: `out`",
    default=None,
    type=str,
)
@click.pass_context
def run(
    ctx,
    benchmark,
    cores,
    dry,
    continue_on_error,
    out_dir,
):
    """Run a benchmark.

    BENCHMARK: Path to benchmark YAML file.

    This command:
    1. Fetches and caches all module repositories
    2. Resolves modules and generates an explicit Snakefile
    3. Runs snakemake on the generated Snakefile

    Examples:

      ob run benchmark.yaml                    # Run full benchmark
      ob run benchmark.yaml --cores 8          # Run with 8 cores
      ob run benchmark.yaml --dry              # Generate Snakefile only
    """
    ctx.ensure_object(dict)

    # Retrieve the global debug flag from the Click context
    debug = ctx.obj.get("DEBUG", False)

    _run_benchmark(
        benchmark_path=benchmark,
        cores=cores,
        dry=dry,
        continue_on_error=continue_on_error,
        out_dir=out_dir,
        debug=debug,
    )


def _run_benchmark(
    benchmark_path,
    cores,
    dry,
    continue_on_error,
    out_dir,
    debug,
):
    """Run a full benchmark."""
    import time

    start_time = time.time()

    # Setup output directory
    out_dir = out_dir if out_dir else "out"
    out_dir_path = Path(out_dir)

    # Load and validate benchmark
    try:
        benchmark_path_abs = Path(benchmark_path).resolve()
        b = BenchmarkExecution(benchmark_path_abs, out_dir_path)
        logger.info("Benchmark YAML file integrity check passed.")
    except PydanticValidationError as e:
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

    # Use clean UI when not in debug mode
    use_clean_ui = not debug

    # Step 1: Populate git cache (fetch all repos)
    _populate_git_cache(b, quiet=use_clean_ui, cores=cores)

    # Step 2: Generate explicit Snakefile
    _generate_explicit_snakefile(
        benchmark=b,
        benchmark_yaml_path=benchmark_path_abs,
        out_dir=out_dir_path,
        debug_mode=False,  # Always generate executable Snakefile
        cores=cores,
        quiet=use_clean_ui,
        start_time=start_time,
    )

    # If dry run, exit here
    if dry:
        software_backend = b.get_benchmark_software_backend()
        hint = ""
        if software_backend == SoftwareBackendEnum.conda:
            hint = " --use-conda"
        logger.info("\nSnakefile generated. To execute:")
        logger.info(f"  cd {out_dir} && snakemake{hint} --cores {cores}")
        sys.exit(0)

    # Step 3: Run snakemake
    _run_snakemake(
        out_dir=out_dir_path,
        cores=cores,
        continue_on_error=continue_on_error,
        software_backend=b.get_benchmark_software_backend(),
        debug=debug,
    )


def _run_snakemake(
    out_dir: Path,
    cores: int,
    continue_on_error: bool,
    software_backend: SoftwareBackendEnum,
    debug: bool,
):
    """
    Run snakemake on the generated Snakefile.

    Args:
        out_dir: Output directory containing the Snakefile
        cores: Number of cores to use
        continue_on_error: Whether to continue on job failures
        software_backend: Software backend (determines --use-conda etc.)
        debug: Whether to enable verbose output
    """
    snakefile_path = out_dir / "Snakefile"

    if not snakefile_path.exists():
        log_error_and_quit(logger, f"Snakefile not found at {snakefile_path}")
        return

    # Build snakemake command
    cmd = ["snakemake", "--snakefile", str(snakefile_path), "--cores", str(cores)]

    # Add software environment flag based on backend
    if software_backend == SoftwareBackendEnum.conda:
        cmd.append("--use-conda")
    elif software_backend == SoftwareBackendEnum.apptainer:
        cmd.append("--use-singularity")
    elif software_backend == SoftwareBackendEnum.docker:
        cmd.append("--use-singularity")
    elif software_backend == SoftwareBackendEnum.envmodules:
        cmd.append("--use-envmodules")
    # host backend: no additional flags needed

    # Add optional flags
    if continue_on_error:
        cmd.append("--keep-going")

    if debug:
        cmd.extend(["--verbose", "--debug"])

    logger.info(f"Running: {' '.join(cmd)}")

    # Change to output directory and run snakemake
    original_dir = os.getcwd()
    try:
        os.chdir(out_dir)
        result = subprocess.run(cmd)

        if result.returncode == 0:
            logger.info("Benchmark run completed successfully.")
            sys.exit(0)
        else:
            logger.error(f"Benchmark run failed with exit code {result.returncode}")
            sys.exit(result.returncode)
    finally:
        os.chdir(original_dir)


def log_error_and_quit(logger, error):
    logger.error(error)
    sys.exit(1)


def _populate_git_cache(
    benchmark: BenchmarkExecution, quiet: bool = False, cores: int = 4
):
    """Populate git cache with repositories from benchmark.

    Args:
        benchmark: BenchmarkExecution instance
        quiet: If True, use progress bar instead of logging
        cores: Number of parallel workers for fetching
    """
    from omnibenchmark.cli.progress import ProgressDisplay

    cache_dir = get_git_cache_dir()

    # Extract unique repositories with their commits
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
        from omnibenchmark.git.cache import parse_repo_url

        repo_path = parse_repo_url(repo_url)
        repo_cache_dir = cache_dir / repo_path

        # Skip fetch if repo exists and we have the commit
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
                repo[commit.encode("ascii")]
                skip_fetch = True
            except (KeyError, Exception):
                pass

        if skip_fetch:
            return (repo_url, True, None)

        if not quiet:
            logger.info(f"Caching {repo_url}")

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


def _generate_explicit_snakefile(
    benchmark: BenchmarkExecution,
    benchmark_yaml_path: Path,
    out_dir: Path,
    nesting_strategy: str = "nested",
    debug_mode: bool = False,
    cores: int = 4,
    quiet: bool = False,
    start_time: float = None,
):
    """
    Generate explicit Snakefile.

    This resolves all modules, dereferences entrypoints, and generates
    a human-readable Snakefile with shell: directives.

    Args:
        benchmark: BenchmarkExecution instance
        benchmark_yaml_path: Path to original benchmark YAML
        out_dir: Output directory
        nesting_strategy: Path nesting strategy - "nested" (default) or "flat"
        debug_mode: If True, generate echo-only debug Snakefile
        cores: Number of parallel workers for module resolution
        quiet: If True, use clean progress UI
        start_time: Start time for timing
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
    benchmark_dir = benchmark_yaml_path.parent

    resolver = ModuleResolver(
        work_base_dir=work_dir,
        output_dir=out_dir,
        software_backend=benchmark.model.get_software_backend(),
        software_environments=benchmark.model.get_software_environments(),
        benchmark_dir=benchmark_dir,
    )

    # Collect all unique modules to resolve
    unique_modules = {}
    for stage in benchmark.model.stages:
        for module in stage.modules:
            if module.id not in unique_modules:
                unique_modules[module.id] = (module, module.software_environment)

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

        if quiet:
            root_logger = logging.getLogger()
            old_handlers = root_logger.handlers[:]
            old_level = root_logger.level
            root_logger.handlers = []
            root_logger.setLevel(logging.CRITICAL + 1)

            omnibenchmark_logger = logging.getLogger("omnibenchmark")
            old_omni_level = omnibenchmark_logger.level
            old_omni_propagate = omnibenchmark_logger.propagate
            old_omni_handlers = omnibenchmark_logger.handlers[:]
            omnibenchmark_logger.handlers = []
            omnibenchmark_logger.setLevel(logging.CRITICAL + 1)
            omnibenchmark_logger.propagate = False

        module_warnings = []
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            if quiet:
                warnings.filterwarnings("ignore")

            try:
                resolved = resolver.resolve(
                    module=module,
                    module_id=module_id,
                    software_environment_id=software_env_id,
                    dirty=False,
                )
                for warning in w:
                    module_warnings.append((module_id, str(warning.message)))

                return (module_id, resolved, None, module_warnings)
            except Exception as e:
                return (module_id, None, str(e), module_warnings)
            finally:
                if quiet:
                    root_logger.handlers = old_handlers
                    root_logger.setLevel(old_level)
                    omnibenchmark_logger.handlers = old_omni_handlers
                    omnibenchmark_logger.setLevel(old_omni_level)
                    omnibenchmark_logger.propagate = old_omni_propagate

    if quiet:
        progress.start_task("Resolving modules", total=len(unique_modules))

    with ThreadPoolExecutor(max_workers=cores) as executor:
        futures = {
            executor.submit(resolve_module_task, mod_id, mod, env_id): mod_id
            for mod_id, (mod, env_id) in unique_modules.items()
        }

        for future in as_completed(futures):
            module_id, resolved, error, module_warnings = future.result()

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

    resolved_nodes = []
    output_to_nodes = {}
    previous_stage_nodes = []

    for stage in benchmark.model.stages:
        current_stage_nodes = []

        for module in stage.modules:
            module_id = module.id

            try:
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
                if stage.inputs and previous_stage_nodes:
                    from itertools import product

                    input_node_combinations = previous_stage_nodes
                    node_combinations = list(
                        product(input_node_combinations, params_list)
                    )
                else:
                    node_combinations = [(None, params) for params in params_list]

                # Create a node for each combination
                for input_node, params in node_combinations:
                    # Check if this combination should be excluded
                    if input_node and module.exclude:
                        if input_node.module_id in module.exclude:
                            logger.debug(
                                f"      Excluding combination: {input_node.module_id} → {module_id}"
                            )
                            continue

                    param_id = f".{params.hash_short()}" if params else ".default"

                    # Build node ID
                    if input_node:
                        node_id = f"{input_node.id}-{stage.id}-{module_id}{param_id}"
                    else:
                        node_id = f"{stage.id}-{module_id}{param_id}"

                    # Resolve inputs from the input_node
                    inputs = {}
                    input_name_mapping = {}
                    base_path = None

                    if input_node:
                        output_id_to_path = {}

                        def get_output_ids_for_node(node):
                            result = {}
                            for s in benchmark.model.stages:
                                if s.id == node.stage_id:
                                    for output_spec in s.outputs:
                                        for output_path in node.outputs:
                                            if s.outputs.index(output_spec) < len(
                                                node.outputs
                                            ):
                                                result[output_spec.id] = node.outputs[
                                                    s.outputs.index(output_spec)
                                                ]
                            return result

                        output_id_to_path.update(get_output_ids_for_node(input_node))

                        def get_ancestor_nodes(node):
                            ancestors = []
                            for prev_node in resolved_nodes:
                                if node.id.startswith(prev_node.id + "-"):
                                    ancestors.append(prev_node)
                                    ancestors.extend(get_ancestor_nodes(prev_node))
                            return ancestors

                        for ancestor in get_ancestor_nodes(input_node):
                            output_id_to_path.update(get_output_ids_for_node(ancestor))

                        if stage.inputs:
                            for input_collection in stage.inputs:
                                for input_id in input_collection.entries:
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

                        if inputs:
                            from pathlib import Path as PathLib

                            first_input = next(iter(inputs.values()))
                            base_path = str(PathLib(first_input).parent)

                    # Determine dataset value
                    if not inputs:
                        dataset_value = module_id
                    elif (
                        input_node
                        and hasattr(input_node, "dataset")
                        and input_node.dataset
                    ):
                        dataset_value = input_node.dataset
                    else:
                        dataset_value = None

                    # Build output paths
                    outputs = []
                    for output_spec in stage.outputs:
                        output_path_template = output_spec.path

                        if dataset_value and "{dataset}" in output_path_template:
                            output_path_template = output_path_template.replace(
                                "{dataset}", dataset_value
                            )

                        if nesting_strategy == "nested":
                            if base_path:
                                output_path = f"{base_path}/{stage.id}/{module_id}/{param_id}/{output_path_template}"
                            else:
                                output_path = f"{stage.id}/{module_id}/{param_id}/{output_path_template}"
                        elif nesting_strategy == "flat":
                            output_path = f"{stage.id}/{module_id}/{param_id}/{output_path_template}"
                        else:
                            raise ValueError(
                                f"Unknown nesting strategy: {nesting_strategy}"
                            )

                        outputs.append(output_path)
                        if output_spec.id not in output_to_nodes:
                            output_to_nodes[output_spec.id] = []
                        output_to_nodes[output_spec.id].append((node_id, output_path))

                    # Determine resource requirements
                    node_resources = None
                    if hasattr(module, "resources") and module.resources:
                        node_resources = module.resources
                    elif hasattr(stage, "resources") and stage.resources:
                        node_resources = stage.resources

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

                if logger.level <= 10:
                    traceback.print_exc()

        previous_stage_nodes = current_stage_nodes

    if not quiet:
        logger.info(f"Created {len(resolved_nodes)} nodes")

    # Resolve metric collectors as regular nodes
    collector_nodes = resolve_metric_collectors(
        metric_collectors=benchmark.model.get_metric_collectors(),
        resolved_nodes=resolved_nodes,
        benchmark=benchmark.model,
        resolver=resolver,
        quiet=quiet,
    )

    resolved_nodes.extend(collector_nodes)

    # Generate Snakefile
    snakefile_path = out_dir / "Snakefile"

    generator = SnakemakeGenerator(
        benchmark_name=benchmark.model.get_name(),
        benchmark_version=benchmark.model.get_version(),
        benchmark_author=benchmark.model.get_author(),
    )

    generator.generate_snakefile(
        nodes=resolved_nodes,
        collectors=[],
        output_path=snakefile_path,
        debug_mode=debug_mode,
    )

    # Save metadata
    save_metadata(
        benchmark_yaml_path=benchmark_yaml_path,
        output_dir=out_dir,
        nodes=resolved_nodes,
        collectors=[],
    )

    if quiet:
        elapsed = time.time() - start_time
        progress.success(
            f"Generated {len(resolved_nodes)} rules in {elapsed:.1f}s in {snakefile_path}"
        )
    else:
        logger.info(f"\nSnakefile generation complete: {snakefile_path}")
        logger.info(f"  Benchmark: {benchmark.model.get_name()}")
        logger.info(f"  Modules: {len(resolved_modules_cache)}")
        logger.info(f"  Nodes: {len(resolved_nodes)}")
