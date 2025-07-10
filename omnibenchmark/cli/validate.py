"""CLI commands for validation of benchmarks and modules."""

import click
from pathlib import Path

from omnibenchmark.cli.utils.logging import logger
from omnibenchmark.cli.utils.validation import validate_benchmark
from omnibenchmark.model.repo import get_repo_hash
from omnibenchmark.benchmark.validation_core import (
    validate_module_files,
    ValidationSeverity,
    ValidationException,
)


@click.command(name="check")
@click.option(
    "--benchmark",
    "-b",
    required=True,
    type=click.Path(exists=True),
    help="Path to benchmark yaml file or benchmark id.",
    envvar="OB_BENCHMARK",
)
@click.option(
    "--no-consistency",
    is_flag=True,
    default=False,
    help="Skip module consistency validation",
)
@click.option(
    "--no-metadata",
    is_flag=True,
    default=False,
    help="Skip metadata validation (CITATION.cff, LICENSE, omnibenchmark.yaml)",
)
@click.option(
    "--no-intermediate",
    is_flag=True,
    default=False,
    help="Skip intermediate results validation",
)
@click.option(
    "--warn",
    is_flag=True,
    default=False,
    help="Convert errors to warnings instead of failing (default: strict mode)",
)
@click.option(
    "--format",
    type=click.Choice(["json", "yaml", "text"], case_sensitive=False),
    default="text",
    help="Output format for validation results",
)
@click.pass_context
def check(ctx, benchmark, no_consistency, no_metadata, no_intermediate, warn, format):
    """Check benchmark configuration and validate modules.

    Validates:
    - Benchmark YAML structure
    - Module metadata (CITATION.cff, LICENSE, omnibenchmark.yaml)
    - License consistency between CITATION.cff and LICENSE files

    TODO: Future validation features:
    - Module consistency validation (--no-consistency flag)
    - Intermediate results validation (--no-intermediate flag)
    - Output formats (--format flag)

    By default runs in strict mode (fails on first error).
    Use --warn to convert errors to warnings and continue processing.
    """
    logger.info(f"Checking benchmark: {benchmark}")

    # Check if any unimplemented flags are used
    if no_consistency or no_intermediate or format != "text":
        raise NotImplementedError(
            "Flags --no-consistency, --no-intermediate, and --format are not yet implemented."
        )

    # Validate the benchmark YAML structure first
    benchmark_obj = validate_benchmark(benchmark, "/tmp")
    if benchmark_obj is None:
        logger.error("Benchmark YAML validation failed")
        ctx.exit(1)

    logger.info("Benchmark YAML validation passed")

    # Skip module validation if --no-metadata is specified
    if no_metadata:
        logger.info("Skipping module metadata validation (--no-metadata specified)")
        logger.info("Benchmark check completed successfully")
        return

    # Get all modules from the benchmark
    modules = benchmark_obj.get_converter().get_modules()
    logger.info(f"Found {len(modules)} modules to validate")

    validation_errors = []
    validation_warnings = []
    modules_processed = 0
    modules_with_issues = 0

    for module_id, module in modules.items():
        logger.info(f"Validating module: {module_id}")
        modules_processed += 1

        # Get repository information
        repo_info = benchmark_obj.get_converter().get_module_repository(module)
        repo_url = getattr(repo_info, "url", None) or getattr(
            repo_info, "repository", None
        )
        commit_hash = (
            getattr(repo_info, "commit_hash", None)
            or getattr(repo_info, "commit", None)
            or getattr(repo_info, "version", None)
        )

        if not repo_url or not commit_hash:
            msg = f"Module {module_id}: Missing repository information"
            if warn:
                validation_warnings.append(msg)
                logger.warning(msg)
            else:
                validation_errors.append(msg)
                logger.error(msg)
                ctx.exit(1)
            continue

        # Check for local repository
        repos_base_dir = Path(".snakemake/repos")
        folder_name = get_repo_hash(repo_url, commit_hash)
        local_repo_path = repos_base_dir / folder_name

        if not local_repo_path.exists():
            msg = f"Module {module_id}: Local repository not found at {local_repo_path}. Run the benchmark first to clone repositories."
            if warn:
                validation_warnings.append(msg)
                logger.warning(msg)
            else:
                validation_errors.append(msg)
                logger.error(msg)
                ctx.exit(1)
            continue

        # Read file contents
        citation_content = None
        license_content = None
        omnibenchmark_content = None
        files_present = {}

        citation_file = local_repo_path / "CITATION.cff"
        license_file = local_repo_path / "LICENSE"
        omnibenchmark_file = local_repo_path / "omnibenchmark.yaml"

        files_present["CITATION.cff"] = citation_file.exists()
        files_present["LICENSE"] = license_file.exists()
        files_present["omnibenchmark.yaml"] = omnibenchmark_file.exists()

        # Read file contents if they exist
        if files_present["CITATION.cff"]:
            try:
                citation_content = citation_file.read_text(encoding="utf-8")
            except Exception as e:
                logger.warning(f"Failed to read CITATION.cff for {module_id}: {e}")

        if files_present["LICENSE"]:
            try:
                license_content = license_file.read_text(encoding="utf-8")
            except Exception as e:
                logger.warning(f"Failed to read LICENSE for {module_id}: {e}")

        if files_present["omnibenchmark.yaml"]:
            try:
                omnibenchmark_content = omnibenchmark_file.read_text(encoding="utf-8")
            except Exception as e:
                logger.warning(
                    f"Failed to read omnibenchmark.yaml for {module_id}: {e}"
                )

        # Validate the module using existing validation functions
        try:
            result = validate_module_files(
                module_id=module_id,
                citation_content=citation_content,
                license_content=license_content,
                omnibenchmark_content=omnibenchmark_content,
                files_present=files_present,
                warn_mode=warn,
                base_path=str(local_repo_path),
            )

            # Process validation results
            if result.errors:
                modules_with_issues += 1
                for error in result.errors:
                    msg = f"Module {module_id}: {error.message}"
                    if warn:
                        validation_warnings.append(msg)
                        logger.warning(msg)
                    else:
                        validation_errors.append(msg)
                        logger.error(msg)
                        if not warn:
                            ctx.exit(1)

            if result.warnings:
                modules_with_issues += 1
                for warning in result.warnings:
                    msg = f"Module {module_id}: {warning.message}"
                    validation_warnings.append(msg)
                    logger.warning(msg)

            if not result.errors and not result.warnings:
                logger.info(f"Module {module_id}: All validations passed ✅")

        except ValidationException as e:
            modules_with_issues += 1
            for issue in e.issues:
                msg = f"Module {module_id}: {issue.message}"
                if issue.severity == ValidationSeverity.ERROR:
                    if warn:
                        validation_warnings.append(msg)
                        logger.warning(msg)
                    else:
                        validation_errors.append(msg)
                        logger.error(msg)
                        if not warn:
                            ctx.exit(1)
                else:
                    validation_warnings.append(msg)
                    logger.warning(msg)

    # Summary
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
