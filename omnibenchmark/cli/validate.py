"""CLI commands for validation of benchmarks and modules."""

import re
import warnings
from pathlib import Path
from typing import Optional

import click
import yaml
from pydantic import ValidationError as PydanticValidationError

from omnibenchmark.core.metadata import (
    ValidationException,
    ValidationSeverity,
    validate_module_files,
)
from omnibenchmark.core.repository_utils import (
    RepositoryManager,
    cleanup_temp_repositories,
)
from omnibenchmark.core.validate import ValidationResult, format_validation_results
from omnibenchmark.logging import logger
from omnibenchmark.model import Benchmark as BenchmarkModel
from omnibenchmark.model.validation import BenchmarkParseError, ValidationError


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


@click.group(name="validate")
@click.pass_context
def validate(ctx):
    """Validate benchmarks and modules."""
    ctx.ensure_object(dict)


@validate.command("plan")
@click.argument("benchmark", type=click.Path(exists=True))
@click.pass_context
def validate_plan(ctx, benchmark: str):
    """Validate benchmark YAML plan structure.

    Checks:
    - YAML syntax is valid
    - Required fields are present
    - Data types are correct
    - References between stages/modules are valid
    - Software environment backends are properly configured
    """
    logger.info(f"Validating benchmark plan: {benchmark}")

    if not (benchmark.endswith(".yaml") or benchmark.endswith(".yml")):
        raise click.UsageError(
            "Invalid benchmark input. Please provide a valid YAML file path."
        )

    # Set up custom warning handler for FutureWarnings
    def custom_warning_handler(
        message, category, filename, lineno, file=None, line=None
    ):
        if category is FutureWarning:
            formatted_msg = _format_warning_message(str(message))
            logger.warning(formatted_msg)
        else:
            # For non-FutureWarnings, use default behavior
            warnings.showwarning(message, category, filename, lineno, file, line)

    # Store old warning handler and install custom one
    old_showwarning = warnings.showwarning
    warnings.showwarning = custom_warning_handler

    # Ensure FutureWarnings are always shown
    warnings.simplefilter("always", FutureWarning)

    try:
        # Load and validate as a Benchmark model
        # This validates YAML syntax, required fields, data types, and references
        benchmark_path = Path(benchmark)
        benchmark_model = BenchmarkModel.from_yaml(benchmark_path)

        if benchmark_model is None:
            logger.error("Error: Failed to parse YAML as a valid OmniBenchmark.")
            ctx.exit(1)

        # Validate execution context (software environments, paths) WITHOUT creating output directory
        # We only need the benchmark directory (parent of YAML file) for resolving relative paths
        benchmark_dir = benchmark_path.parent.absolute()
        benchmark_model.validate_execution_context(benchmark_dir)

        logger.info("✅ Benchmark YAML plan validation passed.")

    except FileNotFoundError:
        logger.error("Error: Benchmark YAML file not found.")
        ctx.exit(1)

    except yaml.YAMLError as e:
        logger.error(f"Error: YAML file format error: {e}")
        ctx.exit(1)

    except PydanticValidationError as e:
        # Handle Pydantic validation errors with detailed field information
        logger.error(f"Error: {format_pydantic_errors(e)}")
        ctx.exit(1)

    except BenchmarkParseError as e:
        logger.error(
            f"Error: Failed to parse YAML as a valid OmniBenchmark: {e.message}"
        )
        ctx.exit(1)

    except ValueError as e:
        logger.error(f"Error: Failed to parse YAML as a valid OmniBenchmark: {str(e)}")
        ctx.exit(1)

    except ValidationError as e:
        # Handle validation errors from validate_execution_context
        logger.error(f"Error: {e}")
        ctx.exit(1)

    except Exception as e:
        logger.error(f"Error: An unexpected error occurred: {e}")
        ctx.exit(1)

    finally:
        # Restore original warning handler
        warnings.showwarning = old_showwarning


@validate.command("module")
@click.argument("path", type=click.Path(exists=True), required=False, default=".")
@click.option(
    "--strict",
    is_flag=True,
    default=False,
    help="Treat warnings as errors and fail on any validation issue.",
)
@click.pass_context
def validate_module(ctx, path, strict, format="summary"):
    """Validate module (metadata, try to run).

    Validates a module repository at PATH. Defaults to current directory if not specified.

    Checks:
    - Required files exist (CITATION.cff, LICENSE, omnibenchmark.yaml)
    - CITATION.cff is valid and has required fields
    - License consistency between CITATION.cff and LICENSE file

    By default, shows warnings but doesn't fail. Use --strict to fail on warnings.
    """
    warn = not strict  # warn=True means lenient mode (default)

    path = Path(path).resolve()
    logger.info(f"Validating module at: {path}")
    if not path.is_dir():
        logger.error(f"Path is not a directory: {path}")
        ctx.exit(1)

    module_id = path.name
    validation_errors = []
    validation_warnings = []

    with RepositoryManager(prefix="omnibenchmark_validate") as repo_manager:
        # Get file contents
        file_contents = repo_manager.get_repository_files(path)
        files_present = repo_manager.get_files_present(path)

        # Validate the module
        try:
            core_result = validate_module_files(
                module_id=module_id,
                citation_content=file_contents.get("citation"),
                license_content=file_contents.get("license"),
                omnibenchmark_content=file_contents.get("omnibenchmark"),
                files_present=files_present,
                warn_mode=warn,
                base_path=str(path),
            )

            validation_result = _convert_validation_result(
                module_id=module_id,
                core_result=core_result,
                repo_url=None,
                commit_hash=None,
                local_repo_exists=True,
                file_contents=file_contents,
            )

            if core_result.errors:
                for error in core_result.errors:
                    msg = f"{error.message}"
                    # Errors should always be errors, even in warn mode
                    validation_errors.append(msg)
                    logger.error(msg)

            if core_result.warnings:
                for warning in core_result.warnings:
                    msg = f"{warning.message}"
                    if warn:
                        validation_warnings.append(msg)
                        logger.warning(msg)
                    else:
                        # In strict mode, treat warnings as errors
                        validation_errors.append(msg)
                        logger.error(msg)

        except ValidationException as e:
            validation_result = _convert_validation_result(
                module_id=module_id,
                core_result=None,
                repo_url=None,
                commit_hash=None,
                local_repo_exists=True,
                file_contents=file_contents,
            )

            for issue in e.issues:
                msg = f"{issue.message}"
                if issue.severity == ValidationSeverity.ERROR:
                    validation_result.add_error(issue.message)
                    # Errors should always be errors
                    validation_errors.append(msg)
                    logger.error(msg)
                else:
                    validation_result.add_warning(issue.message)
                    if warn:
                        validation_warnings.append(msg)
                        logger.warning(msg)
                    else:
                        # In strict mode, treat warnings as errors
                        validation_errors.append(msg)
                        logger.error(msg)

    # Output results
    if format == "summary":
        if validation_errors or validation_warnings:
            # Errors and warnings already logged above
            pass
        else:
            logger.info("✅ Module validation passed")

        if validation_errors:
            ctx.exit(1)
    else:
        validation_results = {module_id: validation_result}
        output = format_validation_results(validation_results, format)
        click.echo(output)

        if validation_errors:
            ctx.exit(1)

    cleanup_temp_repositories()


@validate.command("outputs")
@click.argument("benchmark", type=click.Path(exists=True))
@click.option("--out-dir", default="out", help="Output folder to validate.")
@click.option(
    "--validators-dir",
    default=None,
    help="Validators directory (default: '<plan dir>/<validators.directory>').",
)
@click.option("-c", "--cores", default=1, type=int, help="Snakemake cores.")
@click.option(
    "--force", is_flag=True, default=False, help="Re-validate everything (--forceall)."
)
@click.pass_context
def validate_outputs(ctx, benchmark, out_dir, validators_dir, cores, force):
    """Run output validators against a benchmark's produced files.

    Generates a wildcard Snakefile under ``out/.validation/`` and runs the
    validator declared for each produced output (one per stage/output), in the
    environment named by the plan's ``validators.env``. Pass/fail is recorded per
    output (exit code) for the status report to display.
    """
    rc = run_validators(
        benchmark=benchmark,
        out_dir=out_dir,
        validators_dir=validators_dir,
        cores=cores,
        force=force,
    )
    ctx.exit(rc)


def _resolve_env_directive(model, out_dir: Path, benchmark_dir: Path, env_id: str):
    """Return ``(directive_line, snakemake_backend_flag)`` for the validators env.

    Reuses the module resolver + the same directive mapping as the main run, so
    conda/apptainer/envmodules are activated identically. Host backend -> ("", None).
    """
    from omnibenchmark.backend.resolver import ModuleResolver
    from omnibenchmark.model.benchmark import SoftwareBackendEnum

    backend = model.get_software_backend()
    resolver = ModuleResolver(
        work_base_dir=out_dir / ".modules",
        output_dir=out_dir,
        software_backend=backend,
        software_environments=model.get_software_environments(),
        benchmark_dir=benchmark_dir,
    )
    resolved = resolver._resolve_environment(env_id, "validators")

    def _abs(ref: str) -> str:
        # conda yamls / local SIFs are resolved relative to out_dir, but the
        # validators Snakefile lives in out_dir/.validation — snakemake would
        # resolve a relative path against the wrong dir. Make file refs absolute.
        if "://" in ref or Path(ref).is_absolute():
            return ref
        return str((out_dir / ref).resolve())

    directive = ""
    if resolved is not None:
        bt = resolved.backend_type.value
        if bt == "conda":
            directive = f'conda: "{_abs(resolved.reference)}"'
        elif bt in ("apptainer", "docker"):
            directive = f'container: "{_abs(resolved.reference)}"'
        elif bt == "envmodules":
            directive = f'envmodules: "{resolved.reference}"'

    flag = {
        SoftwareBackendEnum.conda: "--use-conda",
        SoftwareBackendEnum.apptainer: "--use-singularity",
        SoftwareBackendEnum.docker: "--use-singularity",
        SoftwareBackendEnum.envmodules: "--use-envmodules",
    }.get(backend)
    return directive, flag


def run_validators(
    benchmark: str,
    out_dir: str = "out",
    validators_dir: Optional[str] = None,
    cores: int = 1,
    force: bool = False,
) -> int:
    """Discover validators, generate the wildcard Snakefile, run it. Returns rc."""
    import os
    import shutil
    import subprocess
    import sys

    from omnibenchmark.core.status.status import prepare_status
    from omnibenchmark.core.validators import (
        VALIDATION_DIR,
        build_targets,
        discover_validators,
        write_validators_snakefile,
    )

    benchmark_path = Path(benchmark)
    model = BenchmarkModel.from_yaml(benchmark_path)
    cfg = model.get_validators_config()
    if cfg is None:
        logger.error("No 'validators:' section in the benchmark plan.")
        return 1
    if cfg.env not in model.get_software_environments():
        logger.error(
            f"validators env '{cfg.env}' is not defined in software_environments."
        )
        return 1

    out_path = Path(out_dir)
    if not out_path.is_dir():
        logger.error(f"Output dir '{out_path}' not found — run the benchmark first.")
        return 1

    benchmark_dir = benchmark_path.parent.absolute()
    vdir = (
        Path(validators_dir)
        if validators_dir
        else benchmark_dir / (cfg.directory or "validators")
    )
    validators = discover_validators(vdir)
    if not validators:
        logger.warning(f"No validators found under {vdir}; nothing to validate.")
        return 0

    # Match validators to the concrete files the benchmark produced. We reuse the
    # status machinery (rather than globbing) because outputs are nested under
    # hidden .{param} dirs that a glob won't traverse.
    _, filedict, _ = prepare_status(
        str(benchmark_path),
        str(out_path),
        return_all=True,
        cache_dir=Path(".snakemake") / "repos",
    )
    produced = [
        (stage_id, f)
        for stage_id in filedict
        for node in filedict[stage_id]
        for f in filedict[stage_id][node]["observed_output_files"]
    ]
    stage_outputs = {
        stage_id: model.get_stage_outputs(stage_id) for stage_id in model.get_stages()
    }
    targets = build_targets(produced, validators, stage_outputs, out_path)
    if not targets:
        logger.warning(
            "No produced outputs matched a validator; nothing to validate. "
            "(Run the benchmark first, or check validators/{stage}/{output_id}/.)"
        )
        return 0

    directive, backend_flag = _resolve_env_directive(
        model, out_path, benchmark_dir, cfg.env
    )
    snakefile = out_path / VALIDATION_DIR / "Snakefile"
    write_validators_snakefile(targets, directive, snakefile)

    bin_dir = Path(sys.executable).parent
    candidate = bin_dir / "snakemake"
    snakemake_bin = (
        str(candidate)
        if candidate.exists()
        else (shutil.which("snakemake") or "snakemake")
    )
    cmd = [
        snakemake_bin,
        "--snakefile",
        f"{VALIDATION_DIR}/Snakefile",
        "--cores",
        str(cores),
    ]
    if backend_flag:
        cmd.append(backend_flag)
    if force:
        cmd.append("--forceall")

    logger.info(f"Validating {len(targets)} output file(s) -> {snakefile}")
    original_dir = os.getcwd()
    os.chdir(out_path)
    try:
        result = subprocess.run(cmd)
    finally:
        os.chdir(original_dir)
    return result.returncode


def _format_warning_message(message: str) -> str:
    """Format warning message in a concise, colored format.

    Args:
        message: The warning message

    Returns:
        Formatted warning string with yellow color
    """
    line_num = None

    # Simplify common warning messages
    if "should be a string" in message:
        # Extract field name, line number, and value from message
        field_match = re.search(
            r"Field '(\w+)' should be a string(?: \(line (\d+)\))?\. Found \w+ value: ([\d.]+)\.",
            message,
        )
        if field_match:
            field_name = field_match.group(1)
            line_num = field_match.group(2)
            value = field_match.group(3)
            message = f"Field '{field_name}' should be quoted in YAML (found: {value})"
    elif "entries' field in inputs is deprecated" in message:
        # Extract line number if present
        line_match = re.search(r"deprecated(?: \(line (\d+)\))?", message)
        if line_match and line_match.group(1):
            line_num = line_match.group(1)
        message = "Use simple list format for inputs instead of 'entries' field"
    elif "values' parameter format" in message:
        # Extract line number if present
        line_match = re.search(r"deprecated(?: \(line (\d+)\))?", message)
        if line_match and line_match.group(1):
            line_num = line_match.group(1)
        message = "Use dict format for parameters instead of 'values' CLI args"

    # Yellow color for warnings (ANSI code)
    yellow = "\033[93m"
    reset = "\033[0m"

    # Add line number if found
    line_info = f" line {line_num}:" if line_num else ""

    return f"{yellow}[WARN]{reset}{line_info} {message}"


def _convert_validation_result(
    module_id: str,
    core_result,
    repo_url: Optional[str] = None,
    commit_hash: Optional[str] = None,
    local_repo_exists: bool = False,
    file_contents: Optional[dict] = None,
) -> ValidationResult:
    """Convert validation_core result to ValidationResult format."""
    result = ValidationResult(module_id)

    # Basic info
    result.repository_url = repo_url
    result.commit_hash = commit_hash
    result.local_repo_exists = local_repo_exists

    # File existence
    if file_contents:
        result.citation_file_exists = file_contents.get("citation") is not None
        result.license_file_exists = file_contents.get("license") is not None
        result.omnibenchmark_yaml_exists = (
            file_contents.get("omnibenchmark") is not None
        )

    # Extract citation information
    if file_contents and file_contents.get("citation"):
        try:
            citation_data = yaml.safe_load(file_contents["citation"])
            if isinstance(citation_data, dict):
                result.citation_data = citation_data
                result.citation_file_valid = True
                result.citation_license = citation_data.get("license")
                result.citation_has_license = bool(citation_data.get("license"))
                result.citation_has_authors = bool(citation_data.get("authors"))
        except Exception:
            result.citation_file_valid = False

    # Convert errors and warnings
    if core_result:
        for error in core_result.errors:
            result.add_error(error.message)
        for warning in core_result.warnings:
            result.add_warning(warning.message)

    return result
