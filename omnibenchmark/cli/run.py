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

from omnibenchmark.backend._metric_collector import resolve_metric_collectors
from omnibenchmark.backend.resolver import ModuleResolver
from omnibenchmark.backend._manifest import write_run_manifest
from omnibenchmark.backend._metadata import save_metadata
from omnibenchmark.backend.snakemake import SnakemakeGenerator
from omnibenchmark.core import BenchmarkExecution
from omnibenchmark.core._paths import (
    truncate_path_filename,
    collect_path_exclusions,
    is_lineage_excluded,
)
from omnibenchmark.model.params import Params
from omnibenchmark.cli.formatting import pretty_print_parse_error
from omnibenchmark.logging import logger
from omnibenchmark.core import populate_git_cache
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
    "--yes",
    "-y",
    "yes_flag",
    help="Deprecated: accepted for backward compatibility, has no effect.",
    is_flag=True,
    default=False,
    hidden=True,
)
@click.option(
    "--use-remote-storage",
    help="Execute and store results remotely using S3 storage configured in the benchmark YAML.",
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
@click.option(
    "--telemetry",
    is_flag=True,
    default=False,
    help="Emit OTLP telemetry as JSON Lines to stdout (disables Rich progress).",
)
@click.option(
    "--telemetry-output",
    type=click.Path(),
    default=None,
    help="Write telemetry to file instead of stdout. Allows Rich progress to remain active. Implies --telemetry.",
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
    yes_flag,
    use_remote_storage,
    module_filter,
    telemetry,
    telemetry_output,
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
      ob run benchmark.yaml --telemetry         # Emit OTLP/JSONL telemetry to stdout
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
        use_remote_storage=use_remote_storage,
        module_filter=module_filter,
        telemetry=telemetry,
        telemetry_output=telemetry_output,
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
    use_remote_storage=False,
    module_filter=None,
    telemetry=None,
    telemetry_output=None,
    snakemake_args=None,
):
    """Run a full benchmark, or a single-module sub-graph when module_filter is set."""
    start_time = time.time()

    out_dir = out_dir if out_dir else "out"
    out_dir_path = Path(out_dir)

    # --telemetry-output implies --telemetry.
    telemetry = telemetry or bool(telemetry_output)

    # Telemetry to stdout disables Rich progress (they would compete for stdout).
    use_telemetry_stdout = telemetry and not telemetry_output

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

    # Suppress logging when streaming telemetry to stdout (would interleave with JSONL).
    if use_telemetry_stdout:
        logging.getLogger("omnibenchmark").setLevel(logging.WARNING)

    # Setup telemetry emitter early so pre-snakemake phases are traced.
    telemetry_emitter = None
    if telemetry:
        from omnibenchmark.telemetry import TelemetryEmitter

        if telemetry_output:
            telemetry_emitter = TelemetryEmitter(output=Path(telemetry_output))
        else:
            telemetry_emitter = TelemetryEmitter()  # stdout

    # Step 1: Populate git cache (fetch all repos)
    resolution_start_ns = int(time.time() * 1_000_000_000)
    populate_git_cache(b, quiet=use_clean_ui, cores=cores)

    # Step 2: Generate explicit Snakefile
    resolved_nodes = _generate_explicit_snakefile(
        benchmark=b,
        benchmark_yaml_path=benchmark_path_abs,
        out_dir=out_dir_path,
        cores=cores,
        quiet=use_clean_ui,
        start_time=start_time,
        dirty=dirty,
        unpinned=unpinned,
        module_filter=module_filter,
    )

    # Initialize telemetry with resolved nodes and emit module-resolution span
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

    # Build extra Snakemake args, adding S3 remote storage flags if requested
    extra_args = list(snakemake_args or [])
    if use_remote_storage:
        from omnibenchmark.remote.storage import (
            get_storage_from_benchmark,
            remote_storage_snakemake_args,
        )

        # Ensure the S3 bucket exists before running
        get_storage_from_benchmark(b)
        storage_opts = remote_storage_snakemake_args(b)
        for key, value in storage_opts.items():
            if isinstance(value, bool):
                if value:
                    extra_args.append(f"--{key}")
            elif value is not None:
                extra_args.extend([f"--{key}", str(value)])

    # Step 3: Run snakemake
    _run_snakemake(
        out_dir=out_dir_path,
        cores=cores,
        continue_on_error=continue_on_error,
        software_backend=b.get_benchmark_software_backend(),
        debug=debug,
        extra_snakemake_args=extra_args,
        telemetry_emitter=telemetry_emitter,
        benchmark=b,
        resolved_nodes=resolved_nodes,
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
    telemetry_emitter=None,
    benchmark: Optional[BenchmarkExecution] = None,
    resolved_nodes: Optional[list] = None,
):
    """Run snakemake on the generated Snakefile.

    Args:
        out_dir: Output directory containing the Snakefile.
        cores: Number of cores to use.
        continue_on_error: Whether to continue on job failures.
        software_backend: Software backend (determines --use-conda etc.).
        debug: Whether to enable verbose output (shows full output instead of progress).
        extra_snakemake_args: Extra args appended to the snakemake command line.
        telemetry_emitter: Pre-initialized TelemetryEmitter (or None).
        benchmark: BenchmarkExecution instance (for telemetry metadata).
        resolved_nodes: List of ResolvedNode objects (for telemetry manifest).
    """
    from omnibenchmark.progress import ProgressDisplay, InteractiveProgress
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

    # Detect whether the emitter is streaming to stdout. If so, suppress Rich UI.
    use_telemetry_stdout = (
        telemetry_emitter is not None
        and getattr(telemetry_emitter, "_file_handle", None) is sys.stdout
    )

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

    # No Rich console when telemetry is competing for stdout.
    summary_console = None if use_telemetry_stdout else ProgressDisplay().console

    try:
        if debug:
            logger.info(f"Running: {' '.join(cmd)}")
            logger.info(f"Log file: {log_file}")

            with open(log_file, "w") as f:
                result = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)

            with open(log_file, "r", encoding="utf-8", errors="replace") as f:
                print(f.read())
        elif use_telemetry_stdout:
            # Telemetry-to-stdout mode: no Rich progress, just emit telemetry events.
            # use_telemetry_stdout is derived from telemetry_emitter, so it's set here.
            assert telemetry_emitter is not None
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

            current_rule: Optional[str] = None
            completed_jobs = 0
            failed_rules: list[str] = []

            # Buffers for telemetry payloads.
            rule_output_buffer: deque = deque(maxlen=100)
            capturing_error = False
            error_lines: list[str] = []

            setup_buffer: list[str] = []
            first_rule_seen = False
            setup_start_ns = int(time.time() * 1_000_000_000)
            telemetry_emitter.emit_phase_started(
                "environment preparation", "environment_setup"
            )

            assert process.stdout is not None
            with open(log_file, "w") as f:
                for line in process.stdout:
                    f.write(line)
                    f.flush()

                    line_stripped = line.strip()

                    if not first_rule_seen:
                        setup_buffer.append(line_stripped)
                    else:
                        rule_output_buffer.append(line_stripped)

                    if capturing_error:
                        error_lines.append(line_stripped)

                    match = rule_start_pattern.search(line) or localrule_pattern.search(
                        line
                    )
                    if match:
                        if not first_rule_seen:
                            setup_end_ns = int(time.time() * 1_000_000_000)
                            telemetry_emitter.emit_setup_span(
                                "\n".join(setup_buffer),
                                setup_start_ns,
                                setup_end_ns,
                            )
                            first_rule_seen = True

                        current_rule = normalize_rule_name(match.group(1))
                        rule_output_buffer.clear()
                        telemetry_emitter.rule_started(current_rule)

                    if (
                        "Finished job" in line
                        or "Finished jobid" in line
                        or "Nothing to be done" in line
                    ):
                        completed_jobs += 1
                        if current_rule:
                            snakemake_output = "\n".join(rule_output_buffer)
                            rule_log = _read_rule_log(Path("."), current_rule)
                            full_output = (
                                f"{snakemake_output}\n--- Command Output ---\n{rule_log}"
                                if rule_log
                                else snakemake_output
                            )
                            telemetry_emitter.rule_completed(
                                current_rule, stdout=full_output
                            )
                        rule_output_buffer.clear()

                    err_match = job_error_pattern.search(line)
                    if err_match:
                        failed_rule = normalize_rule_name(err_match.group(1))
                        if failed_rule not in failed_rules:
                            failed_rules.append(failed_rule)
                        capturing_error = True
                        error_lines = [line_stripped]

                    if capturing_error and (
                        line_stripped == ""
                        or "Shutting down" in line
                        or "Error executing rule" in line
                    ):
                        if failed_rules:
                            snakemake_output = "\n".join(rule_output_buffer)
                            rule_log = _read_rule_log(Path("."), failed_rules[-1])
                            stderr = (
                                f"{snakemake_output}\n--- Command Output ---\n{rule_log}"
                                if rule_log
                                else (snakemake_output or "\n".join(error_lines))
                            )
                            telemetry_emitter.rule_failed(
                                failed_rules[-1],
                                f"Error in rule {failed_rules[-1]}",
                                stderr=stderr,
                            )
                        capturing_error = False
                        error_lines = []

            process.wait()
            result = process

            # Edge case: no rules ever ran — still emit a setup span.
            if not first_rule_seen and setup_buffer:
                setup_end_ns = int(time.time() * 1_000_000_000)
                telemetry_emitter.emit_setup_span(
                    "\n".join(setup_buffer), setup_start_ns, setup_end_ns
                )

            telemetry_emitter.benchmark_completed(
                success=(result.returncode == 0),
                message=(
                    f"Completed {completed_jobs} jobs"
                    if result.returncode == 0
                    else f"Failed with {len(failed_rules)} error(s)"
                ),
            )
        else:
            # summary_console is only None when use_telemetry_stdout is True,
            # which is handled by the elif branches above; pyright can't see that.
            assert summary_console is not None
            summary_console.print(f"[dim]Running: {' '.join(cmd)}[/dim]")
            summary_console.print(f"[dim]Log file: {log_file}[/dim]")
            summary_console.print()

            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
                encoding="utf-8",  # don't depend on locale
                errors="replace",  # torn multibyte chars on the merged pipe
                # (concurrent --cores>1 jobs) must not kill the runner
            )

            rule_start_pattern = re.compile(r"rule ([\w.]+):")
            job_error_pattern = re.compile(r"Error in rule ([\w.]+):")
            localrule_pattern = re.compile(r"localrule ([\w.]+):")

            progress: InteractiveProgress | None = None
            total_jobs = None
            current_rule = None

            # Telemetry buffers (only used when telemetry_emitter is set).
            rule_output_buffer: deque = deque(maxlen=100)
            capturing_error = False
            error_lines: list[str] = []
            setup_buffer: list[str] = []
            first_rule_seen = False
            setup_start_ns = int(time.time() * 1_000_000_000)
            if telemetry_emitter:
                telemetry_emitter.emit_phase_started(
                    "environment preparation", "environment_setup"
                )

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

                    if telemetry_emitter:
                        if not first_rule_seen:
                            setup_buffer.append(line_stripped)
                        else:
                            rule_output_buffer.append(line_stripped)
                        if capturing_error:
                            error_lines.append(line_stripped)

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
                        if telemetry_emitter and not first_rule_seen:
                            setup_end_ns = int(time.time() * 1_000_000_000)
                            telemetry_emitter.emit_setup_span(
                                "\n".join(setup_buffer),
                                setup_start_ns,
                                setup_end_ns,
                            )
                            first_rule_seen = True

                        current_rule = normalize_rule_name(match.group(1))
                        if telemetry_emitter:
                            rule_output_buffer.clear()
                        if progress:
                            progress.update(current_rule=current_rule)
                        if telemetry_emitter:
                            telemetry_emitter.rule_started(current_rule)

                    match = localrule_pattern.search(line)
                    if match:
                        if telemetry_emitter and not first_rule_seen:
                            setup_end_ns = int(time.time() * 1_000_000_000)
                            telemetry_emitter.emit_setup_span(
                                "\n".join(setup_buffer),
                                setup_start_ns,
                                setup_end_ns,
                            )
                            first_rule_seen = True

                        current_rule = normalize_rule_name(match.group(1))
                        if telemetry_emitter:
                            rule_output_buffer.clear()
                        if progress:
                            progress.update(current_rule=current_rule)
                        if telemetry_emitter:
                            telemetry_emitter.rule_started(current_rule)

                    if (
                        "Finished job" in line
                        or "Finished jobid" in line
                        or "Nothing to be done" in line
                    ):
                        if progress:
                            progress.update(advance=1)
                        if telemetry_emitter and current_rule:
                            snakemake_output = "\n".join(rule_output_buffer)
                            rule_log = _read_rule_log(Path("."), current_rule)
                            full_output = (
                                f"{snakemake_output}\n--- Command Output ---\n{rule_log}"
                                if rule_log
                                else snakemake_output
                            )
                            telemetry_emitter.rule_completed(
                                current_rule, stdout=full_output
                            )
                            rule_output_buffer.clear()

                    match = job_error_pattern.search(line)
                    if match:
                        failed_rule = normalize_rule_name(match.group(1))
                        if progress:
                            progress.add_failed_rule(failed_rule)
                        if telemetry_emitter:
                            capturing_error = True
                            error_lines = [line_stripped]

                    if (
                        telemetry_emitter
                        and capturing_error
                        and (
                            line_stripped == ""
                            or "Shutting down" in line
                            or "Error executing rule" in line
                        )
                    ):
                        if error_lines:
                            err_match = job_error_pattern.search(error_lines[0])
                            if err_match:
                                failed_rule_name = normalize_rule_name(
                                    err_match.group(1)
                                )
                                snakemake_output = "\n".join(rule_output_buffer)
                                rule_log = _read_rule_log(Path("."), failed_rule_name)
                                stderr = (
                                    f"{snakemake_output}\n--- Command Output ---\n{rule_log}"
                                    if rule_log
                                    else (snakemake_output or "\n".join(error_lines))
                                )
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

            # Edge case: no rules ran — still emit a setup span.
            if telemetry_emitter and not first_rule_seen and setup_buffer:
                setup_end_ns = int(time.time() * 1_000_000_000)
                telemetry_emitter.emit_setup_span(
                    "\n".join(setup_buffer), setup_start_ns, setup_end_ns
                )

            completed_jobs = progress.completed if progress else 0
            failed_rules = progress.failed_rules if progress else []

            if progress:
                progress.finish()

            if telemetry_emitter:
                telemetry_emitter.benchmark_completed(
                    success=(result.returncode == 0),
                    message=(
                        f"Completed {completed_jobs} jobs"
                        if result.returncode == 0
                        else f"Failed with {len(failed_rules)} error(s)"
                    ),
                )

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
            if not debug and not use_telemetry_stdout:
                assert summary_console is not None
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
    module_name: Optional[str] = None,
    input_node=None,
    params=None,
) -> TemplateContext:
    """Build a TemplateContext for a node during expansion."""
    provides: dict[str, str] = {}
    module_attrs: dict[str, str] = {
        "id": module_id,
        "stage": stage.id,
        "name": module_name or module_id,
    }

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

    # {name} always resolves to the current module's own ID, never inherited
    provides["name"] = module_id

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


def _lineage_module_ids(input_node, nodes_by_id) -> set:
    """Return the module_ids along input_node's full lineage, inclusive.

    Walks the explicit ``parent_id`` chain rather than matching node-ID string
    prefixes, so it is unaffected by module/stage ids that contain the id
    separator. ``nodes_by_id`` maps node id → node.
    """
    lineage = set()
    node = input_node
    while node is not None:
        lineage.add(node.module_id)
        node = nodes_by_id.get(node.parent_id) if node.parent_id else None
    return lineage


def _generate_explicit_snakefile(
    benchmark: BenchmarkExecution,
    benchmark_yaml_path: Path,
    out_dir: Path,
    nesting_strategy: str = "nested",
    cores: int = 4,
    quiet: bool = False,
    start_time: Optional[float] = None,
    dirty: bool = False,
    unpinned: bool = False,
    module_filter: Optional[str] = None,
):
    """Generate explicit Snakefile from resolved modules."""
    from omnibenchmark.progress import ProgressDisplay
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
    nodes_by_id = {}
    output_to_nodes = {}
    previous_stage_nodes = []
    dag_errors: list[tuple[str, str, str]] = []

    # Lineage-wide exclusion rules, shared with execution-path pruning so the
    # two code paths agree (see core._paths.is_lineage_excluded).
    path_exclusions = collect_path_exclusions(benchmark.model)

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
                    # Exclusions are transitive over the full lineage, not just
                    # the immediate predecessor: prune if any exclusion rule has
                    # both endpoints present along this node's lineage.
                    if input_node:
                        lineage = _lineage_module_ids(input_node, nodes_by_id) | {
                            module_id
                        }
                        if is_lineage_excluded(lineage, path_exclusions):
                            logger.debug(
                                f"      Excluding combination: lineage {lineage} "
                                f"violates an exclusion rule"
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
                                        output_index = s.outputs.index(output_spec)
                                        if output_index < len(node.outputs):
                                            result[output_spec.id] = node.outputs[
                                                output_index
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
                        module_name=getattr(module, "name", None),
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

                        output_path = truncate_path_filename(output_path)

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
                        parent_id=input_node.id if input_node else None,
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
                    nodes_by_id[node.id] = node
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
        api_version=benchmark.model.api_version,
    )

    generator.generate_snakefile(
        nodes=resolved_nodes,
        output_path=snakefile_path,
    )

    save_metadata(
        benchmark_yaml_path=benchmark_yaml_path,
        output_dir=out_dir,
        nodes=resolved_nodes,
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
