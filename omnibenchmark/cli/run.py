"""cli commands related to benchmark/module execution and start"""

import logging
import os
import subprocess
import sys
import time
from collections import deque
from pathlib import Path
from typing import Optional

import click
from pydantic import ValidationError as PydanticValidationError

from omnibenchmark.backend.collector_resolution import resolve_metric_collectors
from omnibenchmark.backend.resolver import ModuleResolver
from omnibenchmark.backend.snakemake_gen import (
    SnakemakeGenerator,
    save_metadata,
    write_run_manifest,
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

        # Make missing field errors more explicit
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
    "--dev",
    help="Allow unpinned branch references on remote repos (resolved to HEAD at run time). Use for development only.",
    is_flag=True,
    default=False,
)
@click.option(
    "--telemetry",
    is_flag=True,
    default=False,
    help=(
        "Emit OTLP telemetry to <out_dir>/telemetry.jsonl and stream it to the "
        "gRPC endpoint set in [telemetry] endpoint (default: localhost:4317). "
        "Requires the telemetry extras (pip install omnibenchmark[telemetry])."
    ),
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
        "first upstream input × parameter expansion for each module. "
        "Not supported for gather stages."
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
    dev,
    telemetry,
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
      ob run benchmark.yaml --dev              # Allow branch refs on remote repos
      ob run benchmark.yaml -m M1              # Dev mode: run only module M1
      ob run benchmark.yaml -- --rerun-triggers mtime  # Pass flags to snakemake
      ob run benchmark.yaml -- --forceall      # Force re-run all rules
    """
    ctx.ensure_object(dict)

    # Retrieve the global debug flag from the Click context
    debug = ctx.obj.get("DEBUG", False)

    # Warn if using dirty or dev mode
    if dirty:
        logger.warning(
            "Running in --dirty mode: local paths with uncommitted changes allowed. "
            "Results may not be reproducible."
        )
    if dev:
        logger.warning(
            "--dev mode: branch refs allowed, results may not be reproducible."
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
        dev=dev,
        telemetry=telemetry,
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
    dev=False,
    telemetry=False,
    module_filter=None,
    snakemake_args=None,
):
    """Run a full benchmark, or a single-module sub-graph when module_filter is set."""
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

    # Setup telemetry emitter early so pre-snakemake phases are traced
    telemetry_emitter = None
    telemetry_endpoint = None
    telemetry_path = None
    if telemetry:
        from omnibenchmark.config import get_telemetry_endpoint
        from omnibenchmark.telemetry import TelemetryEmitter

        out_dir_path.mkdir(parents=True, exist_ok=True)
        telemetry_path = out_dir_path / "telemetry.jsonl"
        telemetry_emitter = TelemetryEmitter(output=telemetry_path)
        telemetry_endpoint = get_telemetry_endpoint()

    # If an endpoint is given, stream telemetry.jsonl to it in a background thread
    _relay_thread = None
    _relay_stop = None
    if telemetry_endpoint and telemetry_emitter:
        import threading

        try:
            import importlib.util as _ilu

            _relay_script = (
                Path(__file__).parent.parent.parent / "scripts" / "telemetry-relay.py"
            )
            _spec = _ilu.spec_from_file_location("telemetry_relay", _relay_script)
            _mod = _ilu.module_from_spec(_spec)
            _spec.loader.exec_module(_mod)
            TelemetryRelay = _mod.TelemetryRelay
        except Exception as _e:
            logger.warning(
                f"Could not load telemetry relay (install opentelemetry-proto grpcio): {_e}"
            )
            TelemetryRelay = None

        if TelemetryRelay is not None:
            import json as _json
            import socket as _socket

            _host, _port_str = telemetry_endpoint.rsplit(":", 1)
            try:
                with _socket.create_connection((_host, int(_port_str)), timeout=2):
                    pass
            except OSError:
                log_error_and_quit(
                    logger,
                    f"Telemetry endpoint {telemetry_endpoint} is not reachable. "
                    "Start the dashboard first (scripts/run_dashboard.sh) or disable --telemetry.",
                )
                return

            _relay_stop = threading.Event()

            def _stream_to_endpoint(path, endpoint, stop_event):
                relay = TelemetryRelay(endpoint=endpoint)
                try:
                    with open(path, "r") as fh:
                        while not stop_event.is_set():
                            line = fh.readline()
                            if not line:
                                stop_event.wait(timeout=0.1)
                                continue
                            line = line.strip()
                            if not line:
                                continue
                            try:
                                relay.send(_json.loads(line))
                            except Exception:
                                pass
                    # Drain any remaining lines after stop
                    for line in fh:
                        line = line.strip()
                        if line:
                            try:
                                relay.send(_json.loads(line))
                            except Exception:
                                pass
                finally:
                    relay.close()

            # Wait for the file to exist (emitter creates it lazily)
            _relay_thread = threading.Thread(
                target=_stream_to_endpoint,
                args=(telemetry_path, telemetry_endpoint, _relay_stop),
                daemon=True,
                name="telemetry-relay",
            )
            _relay_thread.start()
            logger.info(f"Streaming telemetry to {telemetry_endpoint}")

    # Step 1: Populate git cache (fetch all repos)
    resolution_start_ns = int(time.time() * 1_000_000_000)

    _populate_git_cache(b, quiet=use_clean_ui, cores=cores)

    # Step 2: Generate explicit Snakefile
    resolved_nodes = _generate_explicit_snakefile(
        benchmark=b,
        benchmark_yaml_path=benchmark_path_abs,
        out_dir=out_dir_path,
        debug_mode=False,  # Always generate executable Snakefile
        cores=cores,
        quiet=use_clean_ui,
        start_time=start_time,
        dirty=dirty,
        dev=dev,
        module_filter=module_filter,
    )

    # Initialize telemetry with resolved nodes and emit resolution span
    if telemetry_emitter and resolved_nodes:
        stages = [{"id": s.id, "name": s.name or s.id} for s in b.model.stages]
        telemetry_emitter.emit_manifest(
            benchmark_name=b.model.get_name(),
            benchmark_version=b.model.get_version(),
            benchmark_author=b.model.get_author(),
            software_backend=str(b.get_benchmark_software_backend().value),
            cores=cores,
            stages=stages,
            resolved_nodes=resolved_nodes,
        )
        resolution_end_ns = int(time.time() * 1_000_000_000)
        telemetry_emitter.emit_phase_span(
            name="setup: module resolution",
            phase="module_resolution",
            setup_type="resolution",
            output=f"Resolved {len(resolved_nodes)} nodes",
            start_time_ns=resolution_start_ns,
            end_time_ns=resolution_end_ns,
        )

    # Write run manifest (run_id correlates with telemetry trace_id when enabled)
    run_id = telemetry_emitter.run_id if telemetry_emitter else None
    write_run_manifest(output_dir=out_dir_path, run_id=run_id)

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
    try:
        _run_snakemake(
            out_dir=out_dir_path,
            cores=cores,
            continue_on_error=continue_on_error,
            software_backend=b.get_benchmark_software_backend(),
            debug=debug,
            telemetry_emitter=telemetry_emitter,
            benchmark=b,
            resolved_nodes=resolved_nodes,
            extra_snakemake_args=snakemake_args or [],
        )
    finally:
        if _relay_stop is not None:
            _relay_stop.set()
        if _relay_thread is not None:
            _relay_thread.join(timeout=5)


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
    telemetry_emitter=None,
    benchmark: "BenchmarkExecution" = None,
    resolved_nodes: list = None,
    extra_snakemake_args: list = None,
):
    """
    Run snakemake on the generated Snakefile.

    Captures stdout/stderr to log files while showing interactive progress.
    Log files are stored in out_dir/.logs/

    In normal mode, displays an interactive progress panel with:
    - Progress bar showing completed/total jobs
    - Current rule being executed (shortened name)
    - Press 'L' to toggle live log view
    - Press 'Q'/ESC/'P' to return to progress view

    Args:
        out_dir: Output directory containing the Snakefile
        cores: Number of cores to use
        continue_on_error: Whether to continue on job failures
        software_backend: Software backend (determines --use-conda etc.)
        debug: Whether to enable verbose output (shows full output instead of progress)
        telemetry_emitter: Pre-initialized TelemetryEmitter (or None)
        benchmark: BenchmarkExecution instance (for telemetry metadata)
        resolved_nodes: List of ResolvedNode objects (for telemetry manifest)
    """
    from omnibenchmark.cli.progress import ProgressDisplay, InteractiveProgress
    from datetime import datetime
    import re

    snakefile_path = out_dir / "Snakefile"

    if not snakefile_path.exists():
        log_error_and_quit(logger, f"Snakefile not found at {snakefile_path}")
        return

    # Create logs directory (use absolute path since we'll chdir later)
    logs_dir = out_dir.resolve() / ".logs"
    logs_dir.mkdir(parents=True, exist_ok=True)

    # Log file with timestamp (absolute path)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = logs_dir / f"snakemake_{timestamp}.log"

    # Build snakemake command
    # Use just "Snakefile" since we chdir to out_dir before running
    cmd = ["snakemake", "--snakefile", "Snakefile", "--cores", str(cores)]

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

    # Append any extra flags passed after -- on the ob run command line
    if extra_snakemake_args:
        cmd.extend(extra_snakemake_args)

    # Patch the manifest with the exact snakemake invocation (command is now fully built)
    manifest_path = out_dir / ".metadata" / "manifest.json"
    if manifest_path.exists():
        try:
            import json as _json

            _manifest = _json.loads(manifest_path.read_text())
            _manifest["snakemake_cmd"] = cmd
            manifest_path.write_text(_json.dumps(_manifest, indent=2) + "\n")
        except Exception:
            pass

    # Change to output directory
    original_dir = os.getcwd()
    os.chdir(out_dir)

    summary_console = ProgressDisplay().console

    try:
        if debug:
            # In debug mode, show full output directly
            logger.info(f"Running: {' '.join(cmd)}")
            logger.info(f"Log file: {log_file}")

            with open(log_file, "w") as f:
                result = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)

            # Also print to console in debug mode
            with open(log_file, "r") as f:
                print(f.read())
        else:
            # In normal mode, show interactive progress with log toggle
            summary_console.print(f"[dim]Running: {' '.join(cmd)}[/dim]")
            summary_console.print(f"[dim]Log file: {log_file}[/dim]")
            summary_console.print()

            # Start snakemake with output capture
            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
            )

            # Regex patterns for snakemake output
            rule_start_pattern = re.compile(r"rule (\w+):")
            job_error_pattern = re.compile(r"Error in rule (\w+):")
            localrule_pattern = re.compile(r"localrule (\w+):")

            # Interactive progress (starts after we know total jobs)
            progress: InteractiveProgress | None = None
            total_jobs = None
            current_rule = None

            # Buffer to capture output per rule for telemetry
            rule_output_buffer: deque = deque(maxlen=100)  # Last 100 lines per rule
            capturing_error = False
            error_lines = []

            # Track setup phase (before first rule)
            setup_buffer = []
            first_rule_seen = False
            setup_start_ns = int(time.time() * 1_000_000_000)
            if telemetry_emitter:
                telemetry_emitter.emit_phase_started(
                    "environment preparation", "environment_setup"
                )
            setup_status = summary_console.status("Setting up...")
            setup_status.start()

            with open(log_file, "w") as f:
                for line in process.stdout:
                    # Write to log file
                    f.write(line)
                    f.flush()

                    # Check for keyboard input (L=log view, Q/ESC=progress view)
                    if progress:
                        progress.check_keyboard()

                    # Parse for progress info
                    line_stripped = line.strip()

                    # Capture setup output before first rule
                    if not first_rule_seen:
                        setup_buffer.append(line_stripped)
                        if line_stripped and setup_status:
                            setup_status.update(f"[dim]{line_stripped}[/dim]")
                    else:
                        rule_output_buffer.append(line_stripped)

                    # Capture error output after "Error in rule"
                    if capturing_error:
                        error_lines.append(line_stripped)

                    # Check for total job count (appears early in snakemake output)
                    if line_stripped.startswith("total") and total_jobs is None:
                        parts = line_stripped.split()
                        if len(parts) >= 2:
                            try:
                                total_jobs = int(parts[1])
                                # Stop setup status, start interactive progress
                                setup_status.stop()
                                setup_status = None
                                progress = InteractiveProgress(
                                    log_file=log_file, tail_lines=25
                                )
                                progress.start("Running benchmark", total=total_jobs)
                            except ValueError:
                                pass

                    # Check for rule start
                    match = rule_start_pattern.search(line)
                    if match:
                        # Emit setup span when first rule starts
                        if not first_rule_seen and telemetry_emitter:
                            setup_end_ns = int(time.time() * 1_000_000_000)
                            setup_output = "\n".join(setup_buffer)
                            telemetry_emitter.emit_setup_span(
                                setup_output, setup_start_ns, setup_end_ns
                            )
                            first_rule_seen = True

                        current_rule = match.group(1)
                        rule_output_buffer.clear()
                        if progress:
                            progress.update(current_rule=current_rule)
                        if telemetry_emitter:
                            telemetry_emitter.rule_started(current_rule)

                    match = localrule_pattern.search(line)
                    if match:
                        # Emit setup span when first rule starts
                        if not first_rule_seen and telemetry_emitter:
                            setup_end_ns = int(time.time() * 1_000_000_000)
                            setup_output = "\n".join(setup_buffer)
                            telemetry_emitter.emit_setup_span(
                                setup_output, setup_start_ns, setup_end_ns
                            )
                            first_rule_seen = True

                        current_rule = match.group(1)
                        rule_output_buffer.clear()
                        if progress:
                            progress.update(current_rule=current_rule)
                        if telemetry_emitter:
                            telemetry_emitter.rule_started(current_rule)

                    # Check for job completion
                    if "Finished job" in line or "Nothing to be done" in line:
                        if progress:
                            progress.update(advance=1)
                        if telemetry_emitter and current_rule:
                            # Combine snakemake prep output (buffer) with command output (log file)
                            snakemake_output = "\n".join(rule_output_buffer)
                            rule_log = _read_rule_log(Path("."), current_rule)
                            if rule_log:
                                full_output = f"{snakemake_output}\n--- Command Output ---\n{rule_log}"
                            else:
                                full_output = snakemake_output
                            telemetry_emitter.rule_completed(
                                current_rule, stdout=full_output
                            )
                        rule_output_buffer.clear()

                    # Check for errors
                    match = job_error_pattern.search(line)
                    if match:
                        failed_rule = match.group(1)
                        if progress:
                            progress.add_failed_rule(failed_rule)
                        capturing_error = True
                        error_lines = [line_stripped]

                    # End of error block (next rule or empty line after error content)
                    if capturing_error and (
                        line_stripped == ""
                        or "Shutting down" in line
                        or "Error executing rule" in line
                    ):
                        if telemetry_emitter and error_lines:
                            # Get the failed rule from the first error line
                            err_match = job_error_pattern.search(error_lines[0])
                            if err_match:
                                failed_rule_name = err_match.group(1)
                                # Combine snakemake output with per-rule log
                                snakemake_output = "\n".join(rule_output_buffer)
                                rule_log = _read_rule_log(Path("."), failed_rule_name)
                                if rule_log:
                                    stderr = f"{snakemake_output}\n--- Command Output ---\n{rule_log}"
                                else:
                                    stderr = snakemake_output or "\n".join(error_lines)
                                telemetry_emitter.rule_failed(
                                    failed_rule_name,
                                    f"Error in rule {failed_rule_name}",
                                    stderr=stderr,
                                )
                        capturing_error = False
                        error_lines = []

            process.wait()
            result = process

            if setup_status:
                setup_status.stop()

            # Emit setup span if no rules were seen (edge case)
            if not first_rule_seen and telemetry_emitter and setup_buffer:
                setup_end_ns = int(time.time() * 1_000_000_000)
                setup_output = "\n".join(setup_buffer)
                telemetry_emitter.emit_setup_span(
                    setup_output, setup_start_ns, setup_end_ns
                )

            # Get stats before finishing
            completed_jobs = progress.completed if progress else 0
            failed_rules = progress.failed_rules if progress else []

            if progress:
                progress.finish()

            # Emit benchmark completion telemetry
            if telemetry_emitter:
                telemetry_emitter.benchmark_completed(
                    success=(result.returncode == 0),
                    message=f"Completed {completed_jobs} jobs"
                    if result.returncode == 0
                    else f"Failed with {len(failed_rules)} error(s)",
                )

            # Show summary
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
                    for rule in failed_rules[:5]:  # Show first 5
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
    """Populate git cache with repositories from benchmark.

    Args:
        benchmark: BenchmarkExecution instance
        quiet: If True, use progress bar instead of logging
        cores: Number of parallel workers for fetching
    """
    from omnibenchmark.cli.progress import ProgressDisplay

    cache_dir = get_git_cache_dir()

    # Extract unique repositories with their commits (skip local paths)
    repos = {}
    for stage in benchmark.model.stages:
        for module in stage.modules:
            if hasattr(module, "repository") and module.repository:
                repo_url = module.repository.url
                commit = module.repository.commit
                if repo_url and not is_local_path(repo_url):
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
                if repo_url and not is_local_path(repo_url):
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


def _substitute_params_in_path(template: str, params) -> str:
    """Substitute {params.key} placeholders in an output path template.

    Supports dot notation: ``{params.k}`` is replaced by the value of
    parameter ``k``.  Unknown ``{params.*}`` references are left as-is
    so that later stages can still detect unresolved templates.

    Args:
        template: Path template string, e.g. ``"{dataset}_k{params.k}.csv"``
        params: A ``Params`` instance (or None).

    Returns:
        The template with all matching ``{params.*}`` placeholders resolved.
    """
    if params is None or "{params." not in template:
        return template

    import re

    def _replace(match):
        key = match.group(1)
        try:
            return str(params[key])
        except KeyError:
            return match.group(0)  # leave unresolved

    return re.sub(r"\{params\.([^}]+)\}", _replace, template)


def _build_template_context(
    stage,
    module_id: str,
    input_node=None,
    params=None,
) -> TemplateContext:
    """Build a TemplateContext for a node during expansion.

    For entrypoints (no *input_node*): provides come from the stage's own
    ``provides`` list.  When the stage provides a label that also exists as
    a parameter key, the parameter value is used instead of the module ID.
    This lets ``{dataset}`` resolve to the actual parameter value (e.g.
    ``sce_full_Zhengmix4eq``) rather than the module ID literal.

    For downstream map nodes: inherits the input node's provides (full
    ancestor chain) and adds the current stage's provides labels.

    For gather nodes: call with ``input_node=None`` and the stage will
    typically not have ``provides`` — returns a context with only
    ``module_attrs``.
    """
    provides: dict[str, str] = {}
    module_attrs: dict[str, str] = {"id": module_id, "stage": stage.id}

    if input_node is not None:
        # Downstream: inherit ancestor provides (includes "dataset")
        if input_node.template_context is not None:
            provides.update(input_node.template_context.provides)

        # Add own stage's provides labels
        if stage.provides:
            for label in stage.provides:
                if params is not None and label in params:
                    provides[label] = str(params[label])
                else:
                    provides[label] = module_id

        module_attrs["parent.id"] = input_node.module_id
        module_attrs["parent.stage"] = input_node.stage_id
    else:
        # Entrypoint (or gather with no input_node)
        if stage.provides:
            for label in stage.provides:
                # If the providing stage has a parameter with the same name
                # as the label, use the parameter value. This allows
                # path templates like "{dataset}.rds" to resolve to the
                # actual parameter value rather than the module ID.
                if params is not None and label in params:
                    provides[label] = str(params[label])
                else:
                    provides[label] = module_id
                # (output_id is used by gather consumers, not needed here)

        # Entrypoints define the dataset identity
        if not stage.is_gather_stage():
            if params is not None and "dataset" in params:
                provides.setdefault("dataset", str(params["dataset"]))
            else:
                provides.setdefault("dataset", module_id)

    return TemplateContext(provides=provides, module_attrs=module_attrs)


def _build_gather_inputs(
    gather_labels: list[str],
    resolved_nodes: list,
    benchmark_model,
    output_to_nodes: dict[str, list],
) -> tuple[dict[str, str], dict[str, str]]:
    """Build the inputs dict and name mapping for a gather stage.

    For each gathered label, finds the output_id declared in the provider
    stage's ``provides`` map and collects only those paths from
    ``output_to_nodes``.

    Args:
        gather_labels: Labels being gathered (e.g., ["method"])
        resolved_nodes: All resolved nodes so far (unused, kept for signature compat)
        benchmark_model: Benchmark model (for get_provider_stages)
        output_to_nodes: Registry mapping output_id -> [(node_id, path), ...]

    Returns:
        (inputs_dict, input_name_mapping) where:
        - inputs_dict maps "input_0", "input_1", ... to concrete output paths
        - input_name_mapping maps "input_0", ... to the gathered label name
    """
    inputs = {}
    name_mapping = {}
    idx = 0
    for label in gather_labels:
        for provider_stage in benchmark_model.get_provider_stages(label):
            output_id = provider_stage.provides[label]
            for _node_id, output_path in output_to_nodes.get(output_id, []):
                key = f"input_{idx}"
                inputs[key] = output_path
                name_mapping[key] = label
                idx += 1
    return inputs, name_mapping


class GatherOrderingError(Exception):
    """Raised when gather stages appear before their providers in document order."""

    def __init__(self, errors: list[str]):
        self.errors = errors
        super().__init__("\n".join(errors))


def _validate_gather_ordering(stages) -> None:
    """Validate that all provider stages appear before their gather consumers.

    For the MVP we keep document-order processing and simply check that
    every stage referenced via ``provides`` for a gathered label has already
    been seen when we encounter the gather stage.

    Raises:
        GatherOrderingError: if a provider appears after its gather consumer.
    """
    from omnibenchmark.model.benchmark import GatherInput

    # Build provides registry and stage positions in document order
    seen_providers: dict[str, list[str]] = {}  # label -> [stage_id, ...]
    stage_index: dict[str, int] = {}

    for idx, stage in enumerate(stages):
        stage_index[stage.id] = idx
        if stage.provides:
            for label in stage.provides:
                seen_providers.setdefault(label, []).append(
                    stage.id
                )  # provides is now dict

    # Check that every gather reference is satisfied by earlier stages
    errors = []
    for idx, stage in enumerate(stages):
        if not stage.inputs:
            continue
        for inp in stage.inputs:
            if not isinstance(inp, GatherInput):
                continue
            label = inp.gather
            providers = seen_providers.get(label, [])
            for provider_id in providers:
                if stage_index[provider_id] >= idx:
                    errors.append(
                        f"Stage '{stage.id}' gathers '{label}' but provider "
                        f"stage '{provider_id}' appears after it (or is the same stage). "
                        f"Move '{provider_id}' before '{stage.id}' in the YAML."
                    )

    if errors:
        raise GatherOrderingError(errors)


def _generate_explicit_snakefile(
    benchmark: BenchmarkExecution,
    benchmark_yaml_path: Path,
    out_dir: Path,
    nesting_strategy: str = "nested",
    debug_mode: bool = False,
    cores: int = 4,
    quiet: bool = False,
    start_time: float = None,
    dirty: bool = False,
    dev: bool = False,
    module_filter: Optional[str] = None,
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
        dirty: If True, allow unpinned refs (branches). If False, require pinned commits/tags.
        module_filter: If set, prune DAG to the sub-graph needed for this single
            module (its stage and all upstream stages). Within each stage only
            the first input × parameter expansion is kept (see design/005).
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
                    dirty=dirty,
                    dev=dev,
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

    # Collect resolution errors for clear reporting
    resolution_errors = []

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
                resolution_errors.append((module_id, error))
                if not quiet:
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

    # If any modules failed to resolve, show clear error and exit
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

    # Build nodes from resolved modules
    if not quiet:
        logger.info("\nBuilding execution graph...")

    # Validate gather stage ordering: all providers must appear before their consumers
    try:
        _validate_gather_ordering(benchmark.model.stages)
    except GatherOrderingError as e:
        logger.error("Gather ordering error:")
        for err in e.errors:
            logger.error(f"  {err}")
        sys.exit(1)

    # --- Module-filter pruning (ob run -m <module_id>) ---
    # Determine which stages to expand and whether to restrict to first expansion.
    if module_filter:
        # Find the stage that owns the target module
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
        if target_stage.is_gather_stage():
            logger.error(
                f"ob run -m is not supported for gather stages "
                f"(stage '{target_stage.id}' is a gather stage). "
                "See docs/design/005 for details and future work."
            )
            sys.exit(1)

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

    for stage in stages_to_expand:
        current_stage_nodes = []

        # === Gather stage expansion ===
        if stage.is_gather_stage():
            gather_labels = stage.get_gather_labels()

            # Enumerate provider outputs by the output_id declared in provides
            gather_inputs, gather_input_name_mapping = _build_gather_inputs(
                gather_labels, resolved_nodes, benchmark.model, output_to_nodes
            )

            # One node per (module x params) in the gather stage
            for module in stage.modules:
                module_id = module.id

                try:
                    if module_id not in resolved_modules_cache:
                        logger.warning(
                            f"      Module {module_id} not in cache, skipping"
                        )
                        continue

                    resolved_module = resolved_modules_cache[module_id]

                    if module.parameters:
                        params_list = []
                        for param in module.parameters:
                            params_list.extend(Params.expand_from_parameter(param))
                    else:
                        params_list = [None]

                    for params in params_list:
                        param_id = f".{params.hash_short()}" if params else ".default"
                        node_id = f"gather_{stage.id}-{module_id}{param_id}"

                        # Build template context (gather: no provides, only module attrs)
                        ctx = _build_template_context(
                            stage=stage,
                            module_id=module_id,
                        )

                        # Build output paths (no base_path for gather)
                        outputs = []
                        for output_spec in stage.outputs:
                            resolved_path = ctx.substitute(
                                output_spec.path, params=params
                            )
                            output_path = (
                                f"{stage.id}/{module_id}/{param_id}/{resolved_path}"
                            )
                            outputs.append(output_path)
                            output_to_nodes.setdefault(output_spec.id, []).append(
                                (node_id, output_path)
                            )

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
                            inputs=dict(gather_inputs),
                            outputs=outputs,
                            input_name_mapping=dict(gather_input_name_mapping),
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
                            f"      Created {len([n for n in current_stage_nodes if n.module_id == module_id])} "
                            f"gather node(s) for {module_id} "
                            f"(gathering {len(gather_inputs)} inputs)"
                        )

                except Exception as e:
                    logger.error(f"      Failed to resolve {module_id}: {e}")
                    import traceback

                    if logger.isEnabledFor(logging.DEBUG):
                        traceback.print_exc()

            previous_stage_nodes = current_stage_nodes
            continue

        # === Regular (map) stage expansion ===
        # In module-filter mode, determine which modules to expand for this stage:
        #   - target stage: only the named module
        #   - ancestor stages: only the first module (provides a single upstream path)
        if module_filter:
            is_target_stage = stage.id == target_stage.id
            if is_target_stage:
                modules_to_expand = [m for m in stage.modules if m.id == module_filter]
            else:
                # Keep first resolvable module only
                modules_to_expand = stage.modules[:1]
        else:
            modules_to_expand = stage.modules

        for module in modules_to_expand:
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

                # Module-filter mode: keep only the first input × param combination.
                # This produces a minimal smoke-test path while preserving identical
                # output paths to what a full run would generate for the same combo.
                if module_filter:
                    node_combinations = node_combinations[:1]

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

                            # Use the deepest input path as base to preserve full lineage
                            # This ensures outputs nest under the most specific ancestor
                            # e.g., metrics should nest under methods, not datasets
                            deepest_input = max(
                                inputs.values(), key=lambda p: len(PathLib(p).parts)
                            )
                            base_path = str(PathLib(deepest_input).parent)

                    # Build template context
                    ctx = _build_template_context(
                        stage=stage,
                        module_id=module_id,
                        input_node=input_node,
                        params=params,
                    )

                    # Build output paths
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

            except Exception as e:
                logger.error(f"      Failed to resolve {module_id}: {e}")
                import traceback

                if logger.level <= 10:
                    traceback.print_exc()

        previous_stage_nodes = current_stage_nodes

    if not quiet:
        logger.info(f"Created {len(resolved_nodes)} nodes")

    # Resolve metric collectors as regular nodes (skip in -m module mode)
    if not module_filter:
        try:
            collector_nodes = resolve_metric_collectors(
                metric_collectors=benchmark.model.get_metric_collectors(),
                resolved_nodes=resolved_nodes,
                benchmark=benchmark.model,
                resolver=resolver,
                quiet=quiet,
                dirty=dirty,
                dev=dev,
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

    return resolved_nodes
