"""CLI commands for validation of benchmarks and modules."""

from pathlib import Path

import click
import yaml

from omnibenchmark.benchmark import BenchmarkExecution
from omnibenchmark.cli.utils.logging import logger
from omnibenchmark.benchmark.repository_utils import (
    RepositoryManager,
    cleanup_temp_repositories,
    get_module_repository_info,
    resolve_module_repository,
)
from omnibenchmark.benchmark.metadata import (
    validate_module_files,
    ValidationSeverity,
    ValidationException,
)
from omnibenchmark.benchmark.validate import ValidationResult, format_validation_results


def convert_validation_result(
    module_id: str,
    core_result,
    repo_url: str = None,
    commit_hash: str = None,
    local_repo_exists: bool = False,
    file_contents: dict = None,
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


@click.group(name="validate")
@click.pass_context
def validate(ctx):
    """Validate benchmarks, modules, and outputs."""
    ctx.ensure_object(dict)


@validate.command("plan")
@click.option(
    "--benchmark",
    "-b",
    required=True,
    type=click.Path(exists=True),
    help="Path to benchmark yaml file.",
    envvar="OB_BENCHMARK",
)
@click.pass_context
def validate_plan(ctx, benchmark: str):
    """Validate benchmark YAML plan structure.

    Checks:
    - YAML syntax is valid
    - Required fields are present
    - Data types are correct
    - References between stages/modules are valid
    """
    logger.info(f"Validating benchmark plan: {benchmark}")

    if not (benchmark.endswith(".yaml") or benchmark.endswith(".yml")):
        logger.error(
            "Error: Invalid benchmark input. Please provide a valid YAML file path."
        )
        ctx.exit(1)

    try:
        # First check if YAML is parseable
        with open(benchmark, "r") as file:
            yaml.safe_load(file)

        # Then try to load as a Benchmark object (this validates structure)
        b = BenchmarkExecution(benchmark_yaml=Path(benchmark))

        if b is None:
            logger.error("Error: Failed to parse YAML as a valid OmniBenchmark.")
            ctx.exit(1)

        logger.info("✅ Benchmark YAML plan validation passed.")

    except FileNotFoundError:
        logger.error("Error: Benchmark YAML file not found.")
        ctx.exit(1)

    except yaml.YAMLError as e:
        logger.error(f"Error: YAML file format error: {e}")
        ctx.exit(1)

    except ValueError as e:
        logger.error(f"Error: Failed to parse YAML as a valid OmniBenchmark: {str(e)}")
        ctx.exit(1)

    except Exception as e:
        logger.error(f"Error: An unexpected error occurred: {e}")
        ctx.exit(1)


@validate.command("module")
@click.option(
    "--benchmark",
    "-b",
    required=True,
    type=click.Path(exists=True),
    help="Path to benchmark yaml file.",
    envvar="OB_BENCHMARK",
)
@click.option(
    "--module",
    "-m",
    type=str,
    help="Specific module ID to validate. If not specified, validates all modules.",
)
@click.option(
    "--all",
    "validate_all",
    is_flag=True,
    default=False,
    help="Validate all modules (default behavior, kept for clarity).",
)
@click.option(
    "--warn",
    is_flag=True,
    default=False,
    help="Convert errors to warnings instead of failing (default: strict mode).",
)
@click.option(
    "--format",
    type=click.Choice(["summary", "detailed", "json"], case_sensitive=False),
    default="summary",
    help="Output format for validation results.",
)
@click.option(
    "--out",
    type=click.Path(),
    help="Output file to write results to (default: stdout).",
)
@click.pass_context
def validate_module(ctx, benchmark, module, validate_all, warn, format, out):
    """Validate module metadata (CITATION.cff, LICENSE, omnibenchmark.yaml).

    Checks:
    - Repository accessibility
    - Required files exist (CITATION.cff, LICENSE, omnibenchmark.yaml)
    - CITATION.cff is valid and has required fields
    - License consistency between CITATION.cff and LICENSE file

    By default validates all modules. Use -m to validate a specific module.
    """
    logger.info(f"Validating module metadata from: {benchmark}")

    # Load the benchmark
    try:
        b = BenchmarkExecution(benchmark_yaml=Path(benchmark))
        if b is None:
            logger.error("Failed to load benchmark")
            ctx.exit(1)
    except Exception as e:
        logger.error(f"Failed to load benchmark: {e}")
        ctx.exit(1)

    # Get modules to validate
    modules = b.get_converter().get_modules()

    if module:
        # Validate specific module
        if module not in modules:
            logger.error(f"Module '{module}' not found in benchmark")
            logger.info(f"Available modules: {', '.join(modules.keys())}")
            ctx.exit(1)
        modules = {module: modules[module]}
        logger.info(f"Validating module: {module}")
    else:
        logger.info(f"Validating all {len(modules)} modules")

    validation_errors = []
    validation_warnings = []
    modules_processed = 0
    modules_with_issues = 0
    validation_results = {}

    with RepositoryManager(prefix="omnibenchmark_validate") as repo_manager:
        for module_id, mod in modules.items():
            logger.info(f"Validating module: {module_id}")
            modules_processed += 1

            validation_result = ValidationResult(module_id)

            # Get repository information
            repo_url, commit_hash = get_module_repository_info(b, mod)
            validation_result.repository_url = repo_url
            validation_result.commit_hash = commit_hash

            if not repo_url or not commit_hash:
                msg = f"Module {module_id}: Missing repository information"
                validation_result.add_error("Missing repository information")
                if warn:
                    validation_warnings.append(msg)
                    logger.warning(msg)
                else:
                    validation_errors.append(msg)
                    logger.error(msg)
                    ctx.exit(1)
                validation_results[module_id] = validation_result
                continue

            # Resolve repository path
            local_repo_path = resolve_module_repository(b, mod, module_id, repo_manager)

            if local_repo_path is None:
                msg = f"Module {module_id}: Failed to access repository. Run the benchmark first to clone repositories."
                validation_result.add_error("Failed to access repository")
                if warn:
                    validation_warnings.append(msg)
                    logger.warning(msg)
                else:
                    validation_errors.append(msg)
                    logger.error(msg)
                    ctx.exit(1)
                validation_results[module_id] = validation_result
                continue

            # Get file contents
            file_contents = repo_manager.get_repository_files(local_repo_path)
            files_present = repo_manager.get_files_present(local_repo_path)

            # Validate the module
            try:
                core_result = validate_module_files(
                    module_id=module_id,
                    citation_content=file_contents.get("citation"),
                    license_content=file_contents.get("license"),
                    omnibenchmark_content=file_contents.get("omnibenchmark"),
                    files_present=files_present,
                    warn_mode=warn,
                    base_path=str(local_repo_path),
                )

                validation_result = convert_validation_result(
                    module_id=module_id,
                    core_result=core_result,
                    repo_url=repo_url,
                    commit_hash=commit_hash,
                    local_repo_exists=True,
                    file_contents=file_contents,
                )

                if core_result.errors:
                    modules_with_issues += 1
                    for error in core_result.errors:
                        msg = f"Module {module_id}: {error.message}"
                        if warn:
                            validation_warnings.append(msg)
                            logger.warning(msg)
                        else:
                            validation_errors.append(msg)
                            logger.error(msg)
                            if not warn:
                                ctx.exit(1)

                if core_result.warnings:
                    modules_with_issues += 1
                    for warning in core_result.warnings:
                        msg = f"Module {module_id}: {warning.message}"
                        validation_warnings.append(msg)
                        logger.warning(msg)

                if not core_result.errors and not core_result.warnings:
                    logger.info(f"Module {module_id}: All validations passed ✅")

            except ValidationException as e:
                modules_with_issues += 1
                validation_result = convert_validation_result(
                    module_id=module_id,
                    core_result=None,
                    repo_url=repo_url,
                    commit_hash=commit_hash,
                    local_repo_exists=True,
                    file_contents=file_contents,
                )

                for issue in e.issues:
                    msg = f"Module {module_id}: {issue.message}"
                    if issue.severity == ValidationSeverity.ERROR:
                        validation_result.add_error(issue.message)
                        if warn:
                            validation_warnings.append(msg)
                            logger.warning(msg)
                        else:
                            validation_errors.append(msg)
                            logger.error(msg)
                            if not warn:
                                ctx.exit(1)
                    else:
                        validation_result.add_warning(issue.message)
                        validation_warnings.append(msg)
                        logger.warning(msg)

            validation_results[module_id] = validation_result

    # Output results
    if format == "summary":
        logger.info("\nValidation Summary:")
        logger.info(f"  Modules processed: {modules_processed}")
        logger.info(f"  Modules with issues: {modules_with_issues}")
        logger.info(f"  Validation errors: {len(validation_errors)}")
        logger.info(f"  Validation warnings: {len(validation_warnings)}")

        if validation_errors and not warn:
            logger.error("Validation failed with errors")
            ctx.exit(1)
        elif validation_warnings or validation_errors:
            logger.warning("Validation completed with warnings")
        else:
            logger.info("All validations passed successfully ✅")
    else:
        output = format_validation_results(validation_results, format)

        if out:
            try:
                with open(out, "w", encoding="utf-8") as f:
                    f.write(output)
                logger.info(f"Output written to {out}")
            except Exception as e:
                logger.error(f"Failed to write output file: {e}")
                ctx.exit(1)
        else:
            click.echo(output)

        if validation_errors and not warn:
            ctx.exit(1)

    cleanup_temp_repositories()


@validate.command("output")
@click.option(
    "--benchmark",
    "-b",
    required=True,
    type=click.Path(exists=True),
    help="Path to benchmark yaml file.",
    envvar="OB_BENCHMARK",
)
@click.pass_context
def validate_output(ctx, benchmark: str):
    """Validate benchmark outputs (placeholder for future implementation).

    Future checks will include:
    - Output files exist at expected locations
    - Output files match declared schema
    - Output files are accessible and readable
    """
    logger.warning("Output validation is not yet implemented.")
    logger.info("This command is a placeholder for future functionality.")
