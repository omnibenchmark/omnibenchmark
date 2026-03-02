"""cli commands related to benchmark/module execution and start"""

import logging
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Optional

import click
from pydantic import ValidationError as PydanticValidationError

from omnibenchmark.backend.collector_resolution import resolve_metric_collectors
from omnibenchmark.backend.resolver import ModuleResolver
from omnibenchmark.backend.manifest import write_run_manifest
from omnibenchmark.backend.snakemake_gen import (
    SnakemakeGenerator,
    save_metadata,
)
from omnibenchmark.benchmark import BenchmarkExecution
from omnibenchmark.benchmark.params import Params
from omnibenchmark.cli.error_formatting import pretty_print_parse_error
from omnibenchmark.cli.utils.logging import logger
from omnibenchmark.config import get_git_cache_dir
from omnibenchmark.git import get_or_update_cached_repo, is_local_path
from omnibenchmark.model import SoftwareBackendEnum
from omnibenchmark.model.resolved import ResolvedNode, TemplateContext
from omnibenchmark.model.validation import BenchmarkParseError


def format_pydantic_errors(e: PydanticValidationError) -> str:
    """Format Pydantic validation errors to show which fields are missing or invalid."""
    error_lines = ["Validation failed:"]
    for error in e.errors():
        field = " -> ".join(str(loc) for loc in error["loc"])
        msg = error["msg"]
        error_type = error["type"]

        if error_type == "missing":
            error_lines.append(f"  - Missing required field: '{field}'")
        else:
            error_lines.append(f"  - Field '{field}': {msg}")

    return "\n".join(error_lines)


@click.command(
    name="run",
    context_settings=dict(allow_extra_args=True, ignore_unknown_options=True),
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
@click.option(
    "--dirty",
    help="Allow local path module references with uncommitted changes. Use for development only.",
    is_flag=True,
    default=False,
)
@click.option(
    "--unpinned",
    help="Allow unpinned branch references on remote repos (resolved to HEAD at run time). Use for development only.",
    is_flag=True,
    default=False,
)
@click.option(
    "-m",
    "--module",
    "module_filter",
    default=None,
    type=str,
    help=(
        "Run only the sub-graph needed for a single module (development mode). "
        "Prunes all stages after the target module's stage and keeps only the "
        "first upstream input × parameter expansion for each module."
    ),
)
@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)
@click.pass_context
def run(
    ctx,
    benchmark,
    cores,
    dry,
    continue_on_error,
    out_dir,
    dirty,
    unpinned,
    module_filter,
    snakemake_args,
):
    """Run a benchmark.

    BENCHMARK: Path to benchmark YAML file.

    This command:
    1. Fetches and caches all module repositories
    2. Resolves modules and generates an explicit Snakefile
    3. Runs snakemake on the generated Snakefile

    Any arguments after -- are passed directly to snakemake.

    Examples:

      ob run benchmark.yaml                    # Run full benchmark
      ob run benchmark.yaml --cores 8          # Run with 8 cores
      ob run benchmark.yaml --dry              # Generate Snakefile only
      ob run benchmark.yaml --dirty            # Allow local paths with uncommitted changes
      ob run benchmark.yaml --unpinned         # Allow branch refs on remote repos
      ob run benchmark.yaml -m M1              # Dev mode: run only module M1
      ob run benchmark.yaml -- --rerun-triggers mtime  # Pass flags to snakemake
      ob run benchmark.yaml -- --forceall      # Force re-run all rules
    """
    ctx.ensure_object(dict)

    debug = ctx.obj.get("DEBUG", False)

    if dirty:
        logger.warning(
            "Running in --dirty mode: local paths with uncommitted changes allowed. "
            "Results may not be reproducible."
        )
    if unpinned:
        logger.warning(
            "--unpinned mode: branch refs allowed, results may not be reproducible."
        )

    if module_filter:
        logger.warning(
            f"-m {module_filter}: single execution path, first param expansion only."
        )

    _run_benchmark(
        benchmark_path=benchmark,
        cores=cores,
        dry=dry,
        continue_on_error=continue_on_error,
        out_dir=out_dir,
        debug=debug,
        dirty=dirty,
        unpinned=unpinned,
        module_filter=module_filter,
        snakemake_args=list(snakemake_args),
    )


def _run_benchmark(
    benchmark_path,
    cores,
    dry,
    continue_on_error,
    out_dir,
    debug,
    dirty,
    unpinned=False,
    module_filter=None,
    snakemake_args=None,
):
    """Run a full benchmark, or a single-module sub-graph when module_filter is set."""
    start_time = time.time()

    out_dir = out_dir if out_dir else "out"
    out_dir_path = Path(out_dir)

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

    use_clean_ui = not debug

    # Step 1: Populate git cache (fetch all repos)
    _populate_git_cache(b, quiet=use_clean_ui, cores=cores)

    # Step 2: Generate explicit Snakefile
    _generate_explicit_snakefile(
        benchmark=b,
        benchmark_yaml_path=benchmark_path_abs,
        out_dir=out_dir_path,
        debug_mode=False,
        cores=cores,
        quiet=use_clean_ui,
        start_time=start_time,
        dirty=dirty,
        unpinned=unpinned,
        module_filter=module_filter,
    )

    # Write run manifest
    write_run_manifest(output_dir=out_dir_path)

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
        extra_snakemake_args=snakemake_args or [],
    )


def _read_rule_log(out_dir: Path, rule_name: str) -> Optional[str]:
    """Read the per-rule log file if it exists."""
    log_path = out_dir / ".logs" / f"{rule_name}.log"
    if log_path.exists():
        try:
            return log_path.read_text()
        except Exception:
            return None
    return None


def _run_snakemake(
    out_dir: Path,
    cores: int,
    continue_on_error: bool,
    software_backend: SoftwareBackendEnum,
    debug: bool,
    extra_snakemake_args: Optional[list] = None,
):
    """Run snakemake on the generated Snakefile."""
    from omnibenchmark.cli.progress import ProgressDisplay, InteractiveProgress
    from datetime import datetime
    import re

    def normalize_rule_name(name: str) -> str:
        return name.replace(".", "_")

    snakefile_path = out_dir / "Snakefile"

    if not snakefile_path.exists():
        log_error_and_quit(logger, f"Snakefile not found at {snakefile_path}")
        return

    logs_dir = out_dir.resolve() / ".logs"
    logs_dir.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = logs_dir / f"snakemake_{timestamp}.log"

    # Prefer snakemake co-located with the current Python executable (pixi/conda env)
    import shutil
    import sys as _sys

    _bin_dir = Path(_sys.executable).parent
    _snakemake_candidate = _bin_dir / "snakemake"
    if _snakemake_candidate.exists():
        snakemake_bin = str(_snakemake_candidate)
    else:
        snakemake_bin = shutil.which("snakemake") or "snakemake"

    cmd = [snakemake_bin, "--snakefile", "Snakefile", "--cores", str(cores)]

    if software_backend == SoftwareBackendEnum.conda:
        cmd.append("--use-conda")
    elif software_backend == SoftwareBackendEnum.apptainer:
        cmd.append("--use-singularity")
    elif software_backend == SoftwareBackendEnum.docker:
        cmd.append("--use-singularity")
    elif software_backend == SoftwareBackendEnum.envmodules:
        cmd.append("--use-envmodules")

    if continue_on_error:
        cmd.append("--keep-going")

    if debug:
        cmd.extend(["--verbose", "--debug"])

    if extra_snakemake_args:
        cmd.extend(extra_snakemake_args)

    # Patch the manifest with the exact snakemake invocation
    manifest_path = out_dir / ".metadata" / "manifest.json"
    if manifest_path.exists():
        try:
            import json as _json

            _manifest = _json.loads(manifest_path.read_text())
            _manifest["snakemake_cmd"] = cmd
            manifest_path.write_text(_json.dumps(_manifest, indent=2) + "\n")
        except Exception:
            pass

    original_dir = os.getcwd()
    os.chdir(out_dir)

    summary_console = ProgressDisplay().console

    try:
        if debug:
            logger.info(f"Running: {' '.join(cmd)}")
            logger.info(f"Log file: {log_file}")

            with open(log_file, "w") as f:
                result = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)

            with open(log_file, "r") as f:
                print(f.read())
        else:
            summary_console.print(f"[dim]Running: {' '.join(cmd)}[/dim]")
            summary_console.print(f"[dim]Log file: {log_file}[/dim]")
            summary_console.print()

            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
            )

            rule_start_pattern = re.compile(r"rule ([\w.]+):")
            job_error_pattern = re.compile(r"Error in rule ([\w.]+):")
            localrule_pattern = re.compile(r"localrule ([\w.]+):")

            progress: InteractiveProgress | None = None
            total_jobs = None
            current_rule = None

            setup_status = summary_console.status("Setting up...")
            setup_status.start()

            assert process.stdout is not None
            with open(log_file, "w") as f:
                for line in process.stdout:
                    f.write(line)
                    f.flush()

                    if progress:
                        progress.check_keyboard()

                    line_stripped = line.strip()

                    if line_stripped.startswith("total") and total_jobs is None:
                        parts = line_stripped.split()
                        if len(parts) >= 2:
                            try:
                                total_jobs = int(parts[1])
                                if setup_status is not None:
                                    setup_status.stop()
                                setup_status = None
                                progress = InteractiveProgress(
                                    log_file=log_file, tail_lines=25
                                )
                                progress.start("Running benchmark", total=total_jobs)
                            except ValueError:
                                pass

                    match = rule_start_pattern.search(line)
                    if match:
                        current_rule = normalize_rule_name(match.group(1))
                        if progress:
                            progress.update(current_rule=current_rule)

                    match = localrule_pattern.search(line)
                    if match:
                        current_rule = normalize_rule_name(match.group(1))
                        if progress:
                            progress.update(current_rule=current_rule)

                    if (
                        "Finished job" in line
                        or "Finished jobid" in line
                        or "Nothing to be done" in line
                    ):
                        if progress:
                            progress.update(advance=1)

                    match = job_error_pattern.search(line)
                    if match:
                        failed_rule = normalize_rule_name(match.group(1))
                        if progress:
                            progress.add_failed_rule(failed_rule)

            process.wait()
            result = process

            if setup_status:
                setup_status.stop()

            completed_jobs = progress.completed if progress else 0
            failed_rules = progress.failed_rules if progress else []

            if progress:
                progress.finish()

            if result.returncode == 0:
                summary_console.print(
                    f"[green]Completed {completed_jobs} jobs successfully[/green]"
                )
            else:
                if failed_rules:
                    summary_console.print(
                        f"[red]Failed with {len(failed_rules)} error(s)[/red]"
                    )
                    summary_console.print()
                    summary_console.print("[red]Failed rules:[/red]")
                    for rule in failed_rules[:5]:
                        summary_console.print(f"  [red]✗[/red] {rule}")
                    if len(failed_rules) > 5:
                        summary_console.print(f"  ... and {len(failed_rules) - 5} more")
                else:
                    summary_console.print(
                        "[red]Snakemake failed before any rules ran[/red]"
                    )
                summary_console.print()
                summary_console.print(f"[dim]See full log: {log_file}[/dim]")

        if result.returncode == 0:
            if not debug:
                summary_console.print()
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
    """Populate git cache with repositories from benchmark."""
    from omnibenchmark.cli.progress import ProgressDisplay

    cache_dir = get_git_cache_dir()

    repos = {}
    for stage in benchmark.model.stages:
        for module in stage.modules:
            if hasattr(module, "repository") and module.repository:
                repo_url = module.repository.url
                commit = module.repository.commit
                if repo_url and not is_local_path(repo_url):
                    repos[repo_url] = commit

    if (
        hasattr(benchmark.model, "metric_collectors")
        and benchmark.model.metric_collectors
    ):
        for collector in benchmark.model.metric_collectors:
            if hasattr(collector, "repository") and collector.repository:
                repo_url = collector.repository.url
                commit = collector.repository.commit
                if repo_url and not is_local_path(repo_url):
                    repos[repo_url] = commit

    if not repos:
        if not quiet:
            logger.info("No repositories found in benchmark definition")
        return

    if not quiet:
        logger.info(f"Populating cache with {len(repos)} repositories...")

    progress = ProgressDisplay() if quiet else None
    if progress is not None:
        progress.start_task("Fetching repositories to cache", total=len(repos))

    success_count = 0
    failed = []

    from concurrent.futures import ThreadPoolExecutor, as_completed
    import threading

    lock = threading.Lock()

    def fetch_repo_task(repo_url, commit):
        from omnibenchmark.git.cache import parse_repo_url

        repo_path = parse_repo_url(repo_url)
        repo_cache_dir = cache_dir / repo_path

        skip_fetch = False
        if (
            repo_cache_dir.exists()
            and commit
            and len(commit) == 40
            and all(c in "0123456789abcdef" for c in commit.lower())
        ):
            try:
                from dulwich import porcelain
                from dulwich.repo import Repo as DulwichRepo
                from typing import cast as _cast

                repo = _cast(DulwichRepo, porcelain.open_repo(str(repo_cache_dir)))
                repo[commit.encode("ascii")]
                skip_fetch = True
            except (KeyError, Exception):
                pass

        if skip_fetch:
            return (repo_url, True, None)

        if not quiet:
            logger.info(f"Caching {repo_url}")

        try:
            get_or_update_cached_repo(repo_url, cache_dir)
            return (repo_url, True, None)
        except Exception as e:
            return (repo_url, False, str(e))

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

            if progress is not None:
                progress.update(advance=1)

    if progress is not None:
        progress.finish()
        progress.success(f"Cached {success_count}/{len(repos)} repositories")
    else:
        logger.info(f"Successfully cached {success_count}/{len(repos)} repositories")

    if failed:
        logger.warning(f"Failed to cache {len(failed)} repositories:")
        for repo_url, error in failed:
            logger.warning(f"  {repo_url}: {error}")


def _substitute_params_in_path(template: str, params) -> str:
    """Substitute {params.key} placeholders in an output path template."""
    if params is None or "{params." not in template:
        return template

    import re

    def _replace(match):
        key = match.group(1)
        try:
            return str(params[key])
        except KeyError:
            return match.group(0)

    return re.sub(r"\{params\.([^}]+)\}", _replace, template)


def _build_template_context(
    stage,
    module_id: str,
    input_node=None,
    params=None,
) -> TemplateContext:
    """Build a TemplateContext for a node during expansion."""
    provides: dict[str, str] = {}
    module_attrs: dict[str, str] = {"id": module_id, "stage": stage.id}

    stage_provides = getattr(stage, "provides", None)

    if input_node is not None:
        if input_node.template_context is not None:
            provides.update(input_node.template_context.provides)

        if stage_provides:
            for label in stage_provides:
                if params is not None and label in params:
                    provides[label] = str(params[label])
                else:
                    provides[label] = module_id

        module_attrs["parent.id"] = input_node.module_id
        module_attrs["parent.stage"] = input_node.stage_id
    else:
        if stage_provides:
            for label in stage_provides:
                if params is not None and label in params:
                    provides[label] = str(params[label])
                else:
                    provides[label] = module_id

        if params is not None and "dataset" in params:
            provides.setdefault("dataset", str(params["dataset"]))
        else:
            provides.setdefault("dataset", module_id)

    return TemplateContext(provides=provides, module_attrs=module_attrs)


def _select_input_nodes(
    declared_input_ids: list[str],
    output_to_nodes: dict,
    resolved_nodes: list,
    stage_ids_in_order: list[str],
    previous_stage_nodes: list,
) -> list:
    """Return the node list to use as the cartesian expansion base for a stage."""
    if not declared_input_ids:
        return previous_stage_nodes

    providing_stage_id_to_depth: dict[str, int] = {}
    for input_id in declared_input_ids:
        for node_id, _path in output_to_nodes.get(input_id, []):
            node_obj = next(
                (n for n in resolved_nodes if n.id == node_id),
                None,
            )
            if node_obj is not None:
                depth = (
                    stage_ids_in_order.index(node_obj.stage_id)
                    if node_obj.stage_id in stage_ids_in_order
                    else -1
                )
                existing = providing_stage_id_to_depth.get(node_obj.stage_id, -1)
                if depth > existing:
                    providing_stage_id_to_depth[node_obj.stage_id] = depth

    if not providing_stage_id_to_depth:
        return previous_stage_nodes

    deepest_stage_id = max(
        providing_stage_id_to_depth,
        key=providing_stage_id_to_depth.__getitem__,
    )
    return [n for n in resolved_nodes if n.stage_id == deepest_stage_id]


def _satisfies_requires(requires: dict, input_node) -> bool:
    """Return True if the input_node's lineage satisfies all requires constraints."""
    if not input_node.template_context:
        return False
    for label, required_value in requires.items():
        actual_value = input_node.template_context.provides.get(label)
        if actual_value != required_value:
            return False
    return True


def _generate_explicit_snakefile(
    benchmark: BenchmarkExecution,
    benchmark_yaml_path: Path,
    out_dir: Path,
    nesting_strategy: str = "nested",
    debug_mode: bool = False,
    cores: int = 4,
    quiet: bool = False,
    start_time: Optional[float] = None,
    dirty: bool = False,
    unpinned: bool = False,
    module_filter: Optional[str] = None,
):
    """Generate explicit Snakefile from resolved modules."""
    from omnibenchmark.cli.progress import ProgressDisplay
    import time

    progress = ProgressDisplay()
    if start_time is None:
        start_time = time.time()

    if not quiet:
        logger.info("\nGenerating explicit Snakefile...")

    work_dir = out_dir / ".modules"
    benchmark_dir = benchmark_yaml_path.parent

    resolver = ModuleResolver(
        work_base_dir=work_dir,
        output_dir=out_dir,
        software_backend=benchmark.model.get_software_backend(),
        software_environments=benchmark.model.get_software_environments(),
        benchmark_dir=benchmark_dir,
    )

    unique_modules = {}
    for stage in benchmark.model.stages:
        for module in stage.modules:
            cache_key = (stage.id, module.id)
            if cache_key not in unique_modules:
                unique_modules[cache_key] = (module, module.software_environment)

    if not quiet:
        logger.info(f"\nResolving {len(unique_modules)} modules...")

    from concurrent.futures import ThreadPoolExecutor, as_completed
    import threading

    resolved_modules_cache = {}
    resolution_lock = threading.Lock()
    captured_warnings = []

    def resolve_module_task(cache_key, module, software_env_id):
        import warnings

        stage_id, module_id = cache_key
        display_id = f"{stage_id}/{module_id}"

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
                    module_id=display_id,
                    software_environment_id=software_env_id,
                    dirty=dirty,
                    unpinned=unpinned,
                )
                for warning in w:
                    module_warnings.append((display_id, str(warning.message)))

                return (cache_key, resolved, None, module_warnings)
            except Exception as e:
                return (cache_key, None, str(e), module_warnings)
            finally:
                if quiet:
                    root_logger.handlers = old_handlers
                    root_logger.setLevel(old_level)
                    omnibenchmark_logger.handlers = old_omni_handlers
                    omnibenchmark_logger.setLevel(old_omni_level)
                    omnibenchmark_logger.propagate = old_omni_propagate

    if quiet:
        progress.start_task("Resolving modules", total=len(unique_modules))

    resolution_errors = []

    with ThreadPoolExecutor(max_workers=cores) as executor:
        futures = {
            executor.submit(resolve_module_task, cache_key, mod, env_id): cache_key
            for cache_key, (mod, env_id) in unique_modules.items()
        }

        for future in as_completed(futures):
            cache_key, resolved, error, module_warnings = future.result()
            stage_id, module_id = cache_key

            with resolution_lock:
                captured_warnings.extend(module_warnings)

            if error:
                resolution_errors.append((f"{stage_id}/{module_id}", error))
                if not quiet:
                    logger.error(f"Failed to resolve {stage_id}/{module_id}: {error}")
            else:
                with resolution_lock:
                    resolved_modules_cache[cache_key] = resolved

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

    if resolution_errors:
        if quiet:
            progress.console.print()
            progress.console.print("[bold red]Resolution failed:[/bold red]")
            for mod_id, error in resolution_errors:
                progress.console.print(f"  [red]✗[/red] [bold]{mod_id}:[/bold] {error}")
            progress.console.print()
        else:
            logger.error("Resolution failed:")
            for mod_id, error in resolution_errors:
                logger.error(f"  {mod_id}: {error}")
        sys.exit(1)

    if not quiet:
        logger.info("\nBuilding execution graph...")

    # Module-filter pruning (ob run -m <module_id>)
    if module_filter:
        target_stage_idx = None
        for idx, stage in enumerate(benchmark.model.stages):
            if any(m.id == module_filter for m in stage.modules):
                target_stage_idx = idx
                break

        if target_stage_idx is None:
            logger.error(
                f"Module '{module_filter}' not found in benchmark. "
                f"Available modules: "
                + ", ".join(m.id for s in benchmark.model.stages for m in s.modules)
            )
            sys.exit(1)

        target_stage = benchmark.model.stages[target_stage_idx]
        stages_to_expand = benchmark.model.stages[: target_stage_idx + 1]
        logger.info(
            f"Module mode: expanding {len(stages_to_expand)} stage(s) "
            f"(up to and including '{target_stage.id}'), "
            "first expansion only."
        )
    else:
        stages_to_expand = benchmark.model.stages

    resolved_nodes = []
    output_to_nodes = {}
    previous_stage_nodes = []
    dag_errors: list[tuple[str, str, str]] = []

    for stage in stages_to_expand:
        current_stage_nodes = []

        # Module-filter: which modules to expand per stage
        if module_filter:
            is_target_stage = stage.id == target_stage.id
            if is_target_stage:
                modules_to_expand = [m for m in stage.modules if m.id == module_filter]
            else:
                modules_to_expand = stage.modules[:1]
        else:
            modules_to_expand = stage.modules

        for module in modules_to_expand:
            module_id = module.id
            cache_key = (stage.id, module_id)

            try:
                if cache_key not in resolved_modules_cache:
                    logger.warning(f"      Module {module_id} not in cache, skipping")
                    continue

                resolved_module = resolved_modules_cache[cache_key]

                if module.parameters:
                    params_list = []
                    for param in module.parameters:
                        params_list.extend(Params.expand_from_parameter(param))
                else:
                    params_list = [None]

                if stage.inputs and previous_stage_nodes:
                    from itertools import product

                    declared_input_ids = [
                        entry
                        for input_col in stage.inputs
                        if hasattr(input_col, "entries")
                        for entry in input_col.entries
                    ]

                    stage_ids_in_order = [s.id for s in stages_to_expand]
                    input_node_combinations = _select_input_nodes(
                        declared_input_ids=declared_input_ids,
                        output_to_nodes=output_to_nodes,
                        resolved_nodes=resolved_nodes,
                        stage_ids_in_order=stage_ids_in_order,
                        previous_stage_nodes=previous_stage_nodes,
                    )

                    node_combinations = list(
                        product(input_node_combinations, params_list)
                    )
                else:
                    node_combinations = [(None, params) for params in params_list]

                if module_filter:
                    node_combinations = node_combinations[:1]

                for input_node, params in node_combinations:
                    if input_node and module.exclude:
                        if input_node.module_id in module.exclude:
                            logger.debug(
                                f"      Excluding combination: {input_node.module_id} → {module_id}"
                            )
                            continue

                    # Also check: input module may declare that it excludes this method
                    if input_node:
                        input_excludes = benchmark.model.get_module_excludes(
                            input_node.module_id
                        )
                        if input_excludes and module_id in input_excludes:
                            logger.debug(
                                f"      Excluding combination: {input_node.module_id} excludes {module_id}"
                            )
                            continue

                    if input_node and module.requires:
                        if not _satisfies_requires(module.requires, input_node):
                            logger.debug(
                                f"      Skipping combination: requires not satisfied for {module_id} "
                                f"(upstream context: {input_node.template_context.provides})"
                            )
                            continue

                    param_id = f".{params.hash_short()}" if params else ".default"

                    if input_node:
                        node_id = f"{input_node.id}-{stage.id}-{module_id}{param_id}"
                    else:
                        node_id = f"{stage.id}-{module_id}{param_id}"

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
                                if not hasattr(input_collection, "entries"):
                                    continue
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

                            deepest_input = max(
                                inputs.values(), key=lambda p: len(PathLib(p).parts)
                            )
                            base_path = str(PathLib(deepest_input).parent)

                    ctx = _build_template_context(
                        stage=stage,
                        module_id=module_id,
                        input_node=input_node,
                        params=params,
                    )

                    outputs = []
                    for output_spec in stage.outputs:
                        output_path_template = ctx.substitute(
                            output_spec.path, params=params
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

                    node_resources = None
                    if hasattr(module, "resources") and module.resources:
                        node_resources = module.resources
                    elif hasattr(stage, "resources") and stage.resources:
                        node_resources = stage.resources

                    node = ResolvedNode(
                        id=node_id,
                        stage_id=stage.id,
                        module_id=module_id,
                        param_id=param_id,
                        module=resolved_module,
                        parameters=params,
                        inputs=inputs,
                        outputs=outputs,
                        input_name_mapping=input_name_mapping,
                        benchmark_name=benchmark.model.get_name(),
                        benchmark_version=benchmark.model.get_version(),
                        benchmark_author=benchmark.model.get_author(),
                        resources=node_resources,
                        template_context=ctx,
                    )

                    resolved_nodes.append(node)
                    current_stage_nodes.append(node)

                if not quiet:
                    logger.info(
                        f"      Created {len([n for n in current_stage_nodes if n.module_id == module_id])} nodes for {module_id}"
                    )

            except ValueError as e:
                msg = str(e)
                logger.error(f"      Failed to resolve {module_id}: {msg}")
                dag_errors.append((stage.id, module_id, msg))

            except Exception as e:
                logger.error(f"      Failed to resolve {module_id}: {e}")
                import traceback

                if logger.level <= 10:
                    traceback.print_exc()

        previous_stage_nodes = current_stage_nodes

    if not quiet:
        logger.info(f"Created {len(resolved_nodes)} nodes")

    if dag_errors:
        if quiet:
            progress.console.print()
            progress.console.print("[bold red]DAG construction failed:[/bold red]")
            for stage_id, module_id, msg in dag_errors:
                progress.console.print(
                    f"  [red]✗[/red] [bold]{stage_id}/{module_id}:[/bold] {msg}"
                )
            progress.console.print()
        else:
            logger.error("DAG construction failed:")
            for stage_id, module_id, msg in dag_errors:
                logger.error(f"  {stage_id}/{module_id}: {msg}")
        sys.exit(1)

    # Resolve metric collectors (skip in -m module mode)
    if not module_filter:
        try:
            collector_nodes = resolve_metric_collectors(
                metric_collectors=benchmark.model.get_metric_collectors(),
                resolved_nodes=resolved_nodes,
                benchmark=benchmark.model,
                resolver=resolver,
                quiet=quiet,
                dirty=dirty,
                unpinned=unpinned,
            )
        except RuntimeError as e:
            if quiet:
                progress.console.print()
                progress.console.print(
                    f"[bold red]Collector resolution failed:[/bold red] {e}"
                )
                progress.console.print()
            else:
                logger.error(f"Collector resolution failed: {e}")
            sys.exit(1)

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

    return resolved_nodes
