"""cli commands related to benchmark/module execution and start"""

import logging
import os
import re
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
    truncate_filename,
    truncate_path_filename,
    collect_path_exclusions,
    is_lineage_excluded,
)
from omnibenchmark.core._lineage import (
    iter_ancestors,
    join_hash,
    satisfies_requires,
    select_input_nodes,
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
@click.option(
    "--with-capability",
    "with_capability",
    multiple=True,
    metavar="NAME",
    help=(
        "Declare a host capability available on this machine (repeatable). "
        "Modules whose requires_capabilities are not all provided are pruned. "
        "Example: --with-capability gpu --with-capability large_mem"
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
    yes_flag,
    use_remote_storage,
    module_filter,
    telemetry,
    telemetry_output,
    with_capability,
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
        available_capabilities=set(with_capability),
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
    available_capabilities=None,
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
        available_capabilities=available_capabilities,
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
    log_path = out_dir / ".logs" / truncate_filename(f"{rule_name}.log")
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


_PARAM_REF_RE = re.compile(r"\{([A-Za-z0-9_-]+)\.params\.([^{}]+)\}")


def _resolve_param_refs(parent_ctx, params):
    """Substitute ``{label.params.key}`` in parameter VALUES from the lineage.

    Resolved against the PARENT's context, which is already fully resolved, so
    a node can never reference its own label and no fixpoint is needed. Chained
    references (C → B → A) work transitively for free, since what lands on each
    node is the resolved `Params`.

    Only *matched* references are touched — unlike `TemplateContext.substitute`,
    a leftover brace here is a legal literal, not an error. Returns a new
    `Params` when anything changed, otherwise `params` itself: `params_list` is
    built once per module and shared across every input bundle, so it must
    never be mutated in place.
    """
    if params is None:
        return params

    # An entrypoint has no lineage to read, but a reference there must still
    # fail loudly rather than reach the shell verbatim — an empty context
    # raises for every label.
    ctx = parent_ctx if parent_ctx is not None else TemplateContext()
    resolved = None

    for key, value in params.items():
        if not isinstance(value, str) or "{" not in value:
            continue
        whole = _PARAM_REF_RE.fullmatch(value)
        if whole:
            # Whole-value reference keeps the ancestor's native type.
            new_value = ctx.lookup_param(*whole.groups())
        elif _PARAM_REF_RE.search(value):
            new_value = _PARAM_REF_RE.sub(
                lambda m: str(ctx.lookup_param(*m.groups())), value
            )
        else:
            continue
        if resolved is None:
            resolved = Params(params)
        resolved[key] = new_value

    return resolved if resolved is not None else params


def _build_template_context(
    stage,
    module_id: str,
    module_name: Optional[str] = None,
    input_node=None,
    params=None,
    extra_provides=None,
) -> TemplateContext:
    """Build a TemplateContext for a node during expansion.

    `extra_provides` is an optional dict layered on top (authoritative) — gather
    nodes use it to bind their group-key label (the `group_by` stage → group
    value) so output templates can reference it (design 010 §3.3).
    """
    provides: dict[str, str] = {}
    # {label.params.*}: every label binding above also records the params of the
    # node that bound it, so a downstream parameter value can read them.
    provides_params: dict = {}
    module_attrs: dict[str, str] = {
        "id": module_id,
        "stage": stage.id,
        "name": module_name or module_id,
    }

    stage_provides = getattr(stage, "provides", None)

    if input_node is not None:
        if input_node.template_context is not None:
            provides.update(input_node.template_context.provides)
            provides_params.update(input_node.template_context.provides_params)

        if stage_provides:
            for label in stage_provides:
                if params is not None and label in params:
                    provides[label] = str(params[label])
                else:
                    provides[label] = module_id
                provides_params[label] = params

        module_attrs["parent.id"] = input_node.module_id
        module_attrs["parent.stage"] = input_node.stage_id
    else:
        if stage_provides:
            for label in stage_provides:
                if params is not None and label in params:
                    provides[label] = str(params[label])
                else:
                    provides[label] = module_id
                provides_params[label] = params

        if extra_provides is None:
            # Builtin `dataset` is an entrypoint concept. A gather node (the
            # only caller passing extra_provides) binds exactly its group key;
            # leaking `dataset=<gather module id>` would silently shadow the
            # real dataset downstream (design 010 §3.3).
            if params is not None and "dataset" in params:
                provides.setdefault("dataset", str(params["dataset"]))
            else:
                provides.setdefault("dataset", module_id)
            provides_params.setdefault("dataset", params)

    if extra_provides:
        provides.update(extra_provides)

    # {name} always resolves to the current module's own ID, never inherited
    provides["name"] = module_id
    provides_params["name"] = params

    return TemplateContext(
        provides=provides,
        module_attrs=module_attrs,
        provides_params=provides_params,
    )


def _ancestor_module_at_stage(member_id, group_stage, nodes_by_id):
    """The module id of the `group_stage` node a member descends from (the
    group value). Parents-aware: sees through joins and prior gathers. Returns
    None when the member has no ancestor in that stage, or when the ancestry
    crosses a fan-in that yields several distinct modules of that stage —
    grouping would be ambiguous."""
    node = nodes_by_id.get(member_id)
    if node is None:
        return None
    found = {
        n.module_id
        for n in (node, *iter_ancestors(node, nodes_by_id))
        if n.stage_id == group_stage
    }
    if len(found) == 1:
        return found.pop()
    return None


def _expand_gather_stage(
    stage,
    benchmark,
    resolved_modules_cache: dict,
    output_to_nodes: dict,
    nodes_by_id: dict,
    path_exclusions=None,
    module_filter=None,
    target_stage=None,
    available_capabilities=None,
) -> list:
    """Expand a gather stage into ResolvedNodes (design 010 MVP).

    For each `gather` entry, collect every node producing its `from` output id
    and partition them by the ancestor module of the entry's `group_by` stage —
    structural grouping over the lineage chain, no `provides` labels. One node
    is emitted per (group value × module × parameter). The chain is cut
    (`parent_id=None`); outputs are written under `stage.prefix` and registered
    in `output_to_nodes` so downstream stages consume them normally.

    Cross-cutting contracts honoured like the scatter sibling: `module_filter`
    prunes to a single execution path (`ob run -m`), and an exclusion rule
    pairing this stage's module with a member's lineage drops that member from
    that module's gather (transitive exclude at member level).
    """
    nodes: list = []

    # group value -> [(member id, from id, path)]. Insertion order is producer
    # expansion order — deterministic; entries are grouped by from id because
    # the outer loop runs per spec.
    grouped: dict = {}
    # The group-key label: the single group_by stage. Stage.validate_gather
    # guarantees all entries share one axis, so [0] is authoritative. Used to
    # bind the template variable (design 010 §3.3).
    group_label = stage.gather[0].group_by

    for spec in stage.gather:
        producers = output_to_nodes.get(spec.from_, [])
        if not producers:
            raise ValueError(
                f"Stage '{stage.id}' gathers from output id '{spec.from_}', "
                f"which no stage produces (design 010 §3.2)."
            )
        for member_id, path in producers:
            gval = _ancestor_module_at_stage(member_id, spec.group_by, nodes_by_id)
            if gval is None:
                logger.warning(
                    f"      Gather '{stage.id}': member {member_id} has no "
                    f"unambiguous ancestor in stage '{spec.group_by}'; dropped "
                    f"from grouping."
                )
                continue
            grouped.setdefault(gval, []).append((member_id, spec.from_, path))

    if not grouped:
        logger.warning(f"Stage '{stage.id}' gathered no members into any group.")

    # Module-filter: same single-execution-path contract as the scatter path.
    if module_filter:
        if target_stage is not None and stage.id == target_stage.id:
            modules_to_expand = [m for m in stage.modules if m.id == module_filter]
        else:
            modules_to_expand = stage.modules[:1]
    else:
        modules_to_expand, _ = _select_capable_modules(
            stage.modules, module_filter, available_capabilities
        )

    for module in modules_to_expand:
        module_id = module.id
        cache_key = (stage.id, module_id)
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

        combos = [(gval, params) for gval in grouped for params in params_list]
        if module_filter:
            combos = combos[:1]

        for gval, params in combos:
            members = grouped[gval]
            if path_exclusions:
                kept = []
                for member_id, from_id, path in members:
                    member = nodes_by_id.get(member_id)
                    lineage = {module_id}
                    if member is not None:
                        lineage |= _lineage_module_ids(member, nodes_by_id)
                    if is_lineage_excluded(lineage, path_exclusions):
                        logger.debug(
                            f"      Gather '{stage.id}': dropping member "
                            f"{member_id} for module {module_id} "
                            f"(exclusion rule)"
                        )
                        continue
                    kept.append((member_id, from_id, path))
            else:
                kept = members
            if not kept:
                logger.warning(
                    f"      Gather '{stage.id}': group '{gval}' has no "
                    f"members left for module {module_id}; skipping node."
                )
                continue

            param_id = f".{params.hash_short()}" if params else ".default"
            node_id = f"{stage.id}-{module_id}-{gval}{param_id}"

            # Enumerated keys + name mapping: same shape metric collectors
            # use, so _write_gather_shell emits `--<from_id> p1 p2 …`.
            inputs: dict = {}
            input_name_mapping: dict = {}
            for idx, (_mid, from_id, path) in enumerate(kept):
                inputs[f"input_{idx}"] = path
                input_name_mapping[f"input_{idx}"] = from_id
            member_ids = list(dict.fromkeys(mid for mid, _f, _p in kept))

            ctx = _build_template_context(
                stage=stage,
                module_id=module_id,
                module_name=getattr(module, "name", None),
                input_node=None,
                params=params,
                extra_provides={group_label: gval},
            )

            outputs = []
            for output_spec in stage.outputs:
                # Only the group label is bound: a template referencing any
                # other lineage label raises in substitute() — the plan-time
                # error design 010 §3.3 wants, never a silent empty sub.
                tmpl = ctx.substitute(output_spec.path, params=params)
                output_path = truncate_path_filename(
                    f"{stage.prefix}/{gval}/{stage.id}/{module_id}/{param_id}/{tmpl}"
                )
                outputs.append(output_path)
                output_to_nodes.setdefault(output_spec.id, []).append(
                    (node_id, output_path)
                )

            node_resources = None
            if getattr(module, "resources", None):
                node_resources = module.resources
            elif getattr(stage, "resources", None):
                node_resources = stage.resources

            nodes.append(
                ResolvedNode(
                    id=node_id,
                    stage_id=stage.id,
                    module_id=module_id,
                    param_id=param_id,
                    module=resolved_module,
                    parameters=params,
                    parent_id=None,
                    inputs=inputs,
                    outputs=outputs,
                    input_name_mapping=input_name_mapping,
                    benchmark_name=benchmark.model.get_name(),
                    benchmark_version=benchmark.model.get_version(),
                    benchmark_author=benchmark.model.get_author(),
                    resources=node_resources,
                    template_context=ctx,
                    is_gather=True,
                    gathered_from=member_ids,
                    parents=list(member_ids),
                )
            )

    return nodes


def _get_output_ids_for_node(node, benchmark) -> dict:
    """Map each output id declared by `node`'s stage to that node's resolved
    path (index-matched). Used to resolve a downstream stage's declared inputs."""
    result: dict = {}
    for s in benchmark.model.stages:
        if s.id == node.stage_id:
            for output_index, output_spec in enumerate(s.outputs):
                if output_index < len(node.outputs):
                    result[output_spec.id] = node.outputs[output_index]
    return result


def _expansion_segment(param_id: str, members) -> str:
    """The directory segment identifying one expansion of a (stage, module).

    A linear node needs only the parameter hash — its ancestry is already in
    the path prefix. A fan-in node's prefix carries just its deepest input, so
    two joins sharing that branch would collide (Snakemake: AmbiguousRule);
    append the parent-set digest, as the node id already does.
    """
    if len(members) > 1:
        return f"{param_id}-{join_hash(m.id for m in members)}"
    return param_id


def _expand_scatter_stage(
    stage,
    benchmark,
    resolved_modules_cache: dict,
    resolved_nodes: list,
    nodes_by_id: dict,
    output_to_nodes: dict,
    previous_stage_nodes: list,
    stages_to_expand: list,
    path_exclusions,
    nesting_strategy: str,
    module_filter,
    target_stage,
    dag_errors: list,
    quiet: bool,
    available_capabilities: Optional[set] = None,
) -> list:
    """Expand one ordinary (scatter/chain) stage into ResolvedNodes.

    Appends created nodes to `resolved_nodes`/`nodes_by_id` as it builds (the
    ancestry walk reads them mid-expansion) and registers outputs in
    `output_to_nodes`; returns this stage's nodes. Extracted from the former
    inline expansion loop, plus one addition: inputs come as *bundles*
    (`_select_input_bundles`), so a diamond join (#289, api ≥ 0.7.0) expands to
    one node per bundle instead of dropping a branch. The gather sibling is
    `_expand_gather_stage`.
    """
    from itertools import product

    current_stage_nodes: list = []

    # Module-filter: which modules to expand per stage
    if module_filter:
        is_target_stage = stage.id == target_stage.id
        if is_target_stage:
            modules_to_expand = [m for m in stage.modules if m.id == module_filter]
        else:
            modules_to_expand = stage.modules[:1]
    else:
        modules_to_expand, _ = _select_capable_modules(
            stage.modules, module_filter, available_capabilities
        )

    # Input bundles depend only on stage-level state — compute once, not per
    # module. None (vs []) distinguishes "no inputs declared" from "inputs
    # declared but nothing resolvable".
    input_bundles = None
    if stage.inputs and previous_stage_nodes:
        declared_input_ids = [
            entry
            for input_col in stage.inputs
            if hasattr(input_col, "entries")
            for entry in input_col.entries
        ]
        input_bundles = _select_input_bundles(
            declared_input_ids=declared_input_ids,
            output_to_nodes=output_to_nodes,
            resolved_nodes=resolved_nodes,
            stage_ids_in_order=[s.id for s in stages_to_expand],
            previous_stage_nodes=previous_stage_nodes,
            nodes_by_id=nodes_by_id,
        )

    # Input resolution depends only on the bundle (never on module/params):
    # cache (inputs, name mapping, base_path) per member-id tuple.
    resolution_cache: dict = {}

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

            if input_bundles is not None:
                node_combinations = list(product(input_bundles, params_list))
            else:
                node_combinations = [(None, params) for params in params_list]

            if module_filter:
                node_combinations = node_combinations[:1]

            for input_bundle, params in node_combinations:
                # A bundle is the producer node(s) feeding this node: a
                # 1-tuple for a linear node, several for a diamond join
                # (#289). `input_node` is the primary (deepest) branch, kept
                # for the id spine, requires-check and template context.
                members = input_bundle if input_bundle else ()
                input_node = members[0] if members else None
                is_join = len(members) > 1

                # Exclusions are transitive over the full lineage, not just
                # the immediate predecessor: prune if any exclusion rule has
                # both endpoints present along this node's lineage. For a join
                # the lineage is the union over every branch.
                if members:
                    lineage = {module_id}
                    for member in members:
                        lineage |= _lineage_module_ids(member, nodes_by_id)
                    if is_lineage_excluded(lineage, path_exclusions):
                        logger.debug(
                            f"      Excluding combination: lineage {lineage} "
                            f"violates an exclusion rule"
                        )
                        continue

                if input_node and module.requires:
                    if not satisfies_requires(module.requires, input_node):
                        logger.debug(
                            f"      Skipping combination: requires not satisfied for {module_id} "
                            f"(upstream context: {input_node.template_context.provides})"
                        )
                        continue

                # Before the hash: param_id feeds the node id, the output
                # directory segment and the human-readable symlink, so two
                # nodes whose `k` really differs must not share one.
                # ponytail: a bad reference raises out of the per-module `try`
                # below, skipping this module's remaining combinations. Right
                # for a misconfigured reference; if a config ever needs
                # per-combination tolerance, catch here and `continue`.
                params = _resolve_param_refs(
                    input_node.template_context if input_node else None, params
                )

                param_id = f".{params.hash_short()}" if params else ".default"

                if is_join:
                    # No single prefix chain: id is a readable stem plus a
                    # short hash of the (sorted) parents (design 010 §3.9).
                    node_id = f"{stage.id}-{module_id}-{join_hash(m.id for m in members)}{param_id}"
                elif input_node:
                    node_id = f"{input_node.id}-{stage.id}-{module_id}{param_id}"
                else:
                    node_id = f"{stage.id}-{module_id}{param_id}"

                inputs = {}
                input_name_mapping = {}
                base_path = None

                if input_node:
                    bundle_key = tuple(m.id for m in members)
                    cached = resolution_cache.get(bundle_key)
                    if cached is None:
                        output_id_to_path = {}

                        # Collect resolvable outputs from every branch's
                        # producer node and its ancestors (union across the
                        # join). Farthest ancestors first so the NEAREST
                        # producer of a shared output id wins (design 010
                        # §3.1).
                        for member in members:
                            for ancestor in reversed(
                                list(iter_ancestors(member, nodes_by_id))
                            ):
                                output_id_to_path.update(
                                    _get_output_ids_for_node(ancestor, benchmark)
                                )
                            output_id_to_path.update(
                                _get_output_ids_for_node(member, benchmark)
                            )

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
                        resolution_cache[bundle_key] = (
                            inputs,
                            input_name_mapping,
                            base_path,
                        )
                    else:
                        inputs, input_name_mapping, base_path = cached
                    # Fresh dicts per node — cached ones must stay pristine.
                    inputs = dict(inputs)
                    input_name_mapping = dict(input_name_mapping)

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

                    seg = _expansion_segment(param_id, members)
                    if nesting_strategy == "nested":
                        if base_path:
                            output_path = f"{base_path}/{stage.id}/{module_id}/{seg}/{output_path_template}"
                        else:
                            output_path = (
                                f"{stage.id}/{module_id}/{seg}/{output_path_template}"
                            )
                    elif nesting_strategy == "flat":
                        output_path = (
                            f"{stage.id}/{module_id}/{seg}/{output_path_template}"
                        )
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
                    parents=[m.id for m in members] if is_join else [],
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

    return current_stage_nodes


def _stage_ancestry(node, nodes_by_id) -> dict:
    """Map stage_id → node ids over `node` and its full (parents-aware)
    ancestry. The per-stage sets are the join-compatibility fingerprint."""
    out: dict = {node.stage_id: {node.id}}
    for ancestor in iter_ancestors(node, nodes_by_id):
        out.setdefault(ancestor.stage_id, set()).add(ancestor.id)
    return out


def _lineages_consistent(a: dict, b: dict) -> bool:
    """Two ancestries may be joined iff they agree wherever they overlap: for
    every stage present in both, they share at least one node. Divergence at
    any common stage (different nodes of the same stage) is a cross-lineage
    pairing and must not be joined — sharing a root is not enough (#289)."""
    for stage_id, ids in a.items():
        other = b.get(stage_id)
        if other is not None and ids.isdisjoint(other):
            return False
    return True


def _select_input_bundles(
    declared_input_ids: list,
    output_to_nodes: dict,
    resolved_nodes: list,
    stage_ids_in_order: list,
    previous_stage_nodes: list,
    nodes_by_id: dict,
) -> list:
    """Return input *bundles* — each a tuple of producer nodes to join into one
    downstream node (issue #289).

    Fast path (linear / single producing stage): every bundle is a 1-tuple and
    the result is identical to wrapping `select_input_nodes` — no behaviour
    change for existing benchmarks. Join path (a declared input is produced only
    on a divergent branch): the anchor is paired with every producer of each
    missing input whose ancestry is *consistent* with the anchor's (and with
    the other partners') — they agree at every shared stage, not merely at the
    root — and the cartesian product of those partners yields one bundle per
    emitted node.
    """
    anchors = select_input_nodes(
        declared_input_ids,
        output_to_nodes,
        resolved_nodes,
        stage_ids_in_order,
        previous_stage_nodes,
        nodes_by_id,
    )
    if not declared_input_ids:
        return [(a,) for a in anchors]

    producers: dict = {}
    for input_id in declared_input_ids:
        producers[input_id] = [
            nodes_by_id[node_id]
            for node_id, _path in output_to_nodes.get(input_id, [])
            if node_id in nodes_by_id
        ]

    from itertools import product as _product

    ancestry_cache: dict = {}

    def _ancestry(n):
        cached = ancestry_cache.get(n.id)
        if cached is None:
            cached = _stage_ancestry(n, nodes_by_id)
            ancestry_cache[n.id] = cached
        return cached

    bundles: list = []
    for anchor in anchors:
        anchor_map = _ancestry(anchor)
        anchor_ids = {nid for ids in anchor_map.values() for nid in ids}
        missing = [
            input_id
            for input_id in declared_input_ids
            if not any(p.id in anchor_ids for p in producers.get(input_id, []))
        ]
        if not missing:
            bundles.append((anchor,))  # linear / covered — fast path
            continue
        # Join: for each uncovered input, take lineage-consistent partners.
        partner_lists = []
        for input_id in missing:
            partners = [
                p
                for p in producers[input_id]
                if _lineages_consistent(anchor_map, _ancestry(p))
            ]
            partner_lists.append(partners)
        if not all(partner_lists):
            unmatched = [i for i, pl in zip(missing, partner_lists) if not pl]
            logger.warning(
                f"      No lineage-compatible producer for input(s) {unmatched} "
                f"of {anchor.id}; the input will not resolve."
            )
            bundles.append((anchor,))
            continue
        emitted = 0
        for combo in _product(*partner_lists):
            # Partners must also be consistent with EACH OTHER (3-way joins).
            if all(
                _lineages_consistent(_ancestry(x), _ancestry(y))
                for i, x in enumerate(combo)
                for y in combo[i + 1 :]
            ):
                bundles.append((anchor,) + combo)
                emitted += 1
        if not emitted:
            logger.warning(
                f"      No mutually consistent partner combination for "
                f"{anchor.id}; the inputs will not resolve."
            )
            bundles.append((anchor,))
    return bundles


def _lineage_module_ids(input_node, nodes_by_id) -> set:
    """Return the module_ids along input_node's GATING lineage, inclusive.

    Walks the explicit ``parent_id`` chain plus join ``parents`` edges (so
    exclusions see every branch of a diamond, #289), but does NOT cross a
    gather cut — the partition forgets (design 010 §3.3), and one excluded
    member among hundreds gathered must not poison every downstream node.
    """
    lineage = {input_node.module_id}
    for ancestor in iter_ancestors(input_node, nodes_by_id, through_gather=False):
        lineage.add(ancestor.module_id)
    return lineage


def _module_capabilities_met(module, available_capabilities) -> bool:
    """True unless the module declares a capability the host did not provide.

    Modules with no `requires_capabilities` always pass. `available_capabilities`
    of None is treated as the empty set (no host facts declared).
    """
    required = getattr(module, "requires_capabilities", None)
    if not required:
        return True
    return set(required).issubset(available_capabilities or set())


def _select_capable_modules(modules, module_filter, available_capabilities):
    """Partition modules into (kept, pruned) by the host-capability gate.

    A module is pruned when its `requires_capabilities` are not all provided.
    `-m/--module` (dev mode) bypasses the gate entirely — the user asked for a
    specific module explicitly, so nothing is silently dropped from under them.
    """
    if module_filter:
        return list(modules), []
    kept, pruned = [], []
    for module in modules:
        (
            kept if _module_capabilities_met(module, available_capabilities) else pruned
        ).append(module)
    return kept, pruned


def _capability_prune_summary(pruned_modules, available_capabilities) -> str:
    """One-line warning naming what the capability gate dropped and how to keep it."""
    ids = sorted({m.id for m in pruned_modules})
    missing = sorted(
        {
            cap
            for m in pruned_modules
            for cap in m.requires_capabilities
            if cap not in (available_capabilities or set())
        }
    )
    flags = " ".join(f"--with-capability {cap}" for cap in missing)
    return (
        f"{len(ids)} module(s) pruned by capability gate: {', '.join(ids)}. "
        f"Rerun with {flags} to include them."
    )


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
    available_capabilities: Optional[set] = None,
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

    # Capability gate: drop modules whose required host capabilities are not all
    # provided, before any checkout/env resolution. -m/--module bypasses it.
    unique_modules = {}
    pruned_modules = []
    for stage in benchmark.model.stages:
        kept, pruned = _select_capable_modules(
            stage.modules, module_filter, available_capabilities
        )
        pruned_modules.extend(pruned)
        for module in pruned:
            logger.warning(
                f"Pruning module '{module.id}': requires capabilities "
                f"{module.requires_capabilities}, host provides "
                f"{sorted(available_capabilities or [])}."
            )
        for module in kept:
            cache_key = (stage.id, module.id)
            if cache_key not in unique_modules:
                unique_modules[cache_key] = (module, module.software_environment)

    if pruned_modules:
        logger.warning(
            _capability_prune_summary(pruned_modules, available_capabilities)
        )

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

    # Expand stages in topological (dependency) order rather than declaration
    # order. select_input_nodes can only pick a parent among already-expanded
    # stages, so a stage must be expanded after every stage producing its inputs.
    # Sorting here makes declaration order irrelevant: a plan that declares a
    # stage before an upstream producer still resolves, as long as that producer
    # is genuinely upstream (not on a divergent branch -- see the diamond guard
    # in model/validation.py). The sort is stable for already-ordered plans, so
    # well-formed benchmarks are unaffected. See
    # https://github.com/omnibenchmark/omnibenchmark/issues/289.
    from omnibenchmark.core._graph import build_stage_dag, compute_stage_order

    topo_order = compute_stage_order(build_stage_dag(benchmark.model))
    topo_index = {stage_id: i for i, stage_id in enumerate(topo_order)}
    topo_sorted_stages = sorted(
        benchmark.model.stages,
        key=lambda s: topo_index.get(s.id, len(topo_index)),
    )

    # Module-filter pruning (ob run -m <module_id>)
    if module_filter:
        target_stage = next(
            (
                stage
                for stage in benchmark.model.stages
                if any(m.id == module_filter for m in stage.modules)
            ),
            None,
        )

        if target_stage is None:
            logger.error(
                f"Module '{module_filter}' not found in benchmark. "
                f"Available modules: "
                + ", ".join(m.id for s in benchmark.model.stages for m in s.modules)
            )
            sys.exit(1)

        target_pos = topo_index.get(target_stage.id, len(topo_index))
        stages_to_expand = [
            stage
            for stage in topo_sorted_stages
            if topo_index.get(stage.id, len(topo_index)) <= target_pos
        ]
        logger.info(
            f"Module mode: expanding {len(stages_to_expand)} stage(s) "
            f"(up to and including '{target_stage.id}'), "
            "first expansion only."
        )
    else:
        stages_to_expand = topo_sorted_stages

    resolved_nodes = []
    nodes_by_id = {}
    output_to_nodes = {}
    previous_stage_nodes = []
    dag_errors: list[tuple[str, str, str]] = []

    # Lineage-wide exclusion rules, shared with execution-path pruning so the
    # two code paths agree (see core._paths.is_lineage_excluded).
    path_exclusions = collect_path_exclusions(benchmark.model)

    for stage in stages_to_expand:
        if stage.gather:
            # Gather stages fan in instead of chaining (design 010). Plan-time
            # errors (zero-producer `from`, unbound template label) feed the
            # same dag_errors report as scatter failures, not a raw traceback.
            try:
                current_stage_nodes = _expand_gather_stage(
                    stage=stage,
                    benchmark=benchmark,
                    resolved_modules_cache=resolved_modules_cache,
                    output_to_nodes=output_to_nodes,
                    nodes_by_id=nodes_by_id,
                    path_exclusions=path_exclusions,
                    module_filter=module_filter,
                    target_stage=target_stage if module_filter else None,
                    available_capabilities=available_capabilities,
                )
            except ValueError as e:
                logger.error(f"      Failed to expand gather stage {stage.id}: {e}")
                dag_errors.append((stage.id, "<gather>", str(e)))
                current_stage_nodes = []
            for node in current_stage_nodes:
                resolved_nodes.append(node)
                nodes_by_id[node.id] = node
        else:
            current_stage_nodes = _expand_scatter_stage(
                stage=stage,
                benchmark=benchmark,
                resolved_modules_cache=resolved_modules_cache,
                resolved_nodes=resolved_nodes,
                nodes_by_id=nodes_by_id,
                output_to_nodes=output_to_nodes,
                previous_stage_nodes=previous_stage_nodes,
                stages_to_expand=stages_to_expand,
                path_exclusions=path_exclusions,
                nesting_strategy=nesting_strategy,
                module_filter=module_filter,
                target_stage=target_stage if module_filter else None,
                dag_errors=dag_errors,
                quiet=quiet,
                available_capabilities=available_capabilities,
            )

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
