"""CLI commands for validation of benchmarks and modules."""

import click

from omnibenchmark.cli.utils.logging import logger
from omnibenchmark.cli.utils.validation import validate_benchmark
from omnibenchmark.benchmark.repository_utils import (
    RepositoryManager,
    cleanup_temp_repositories,
    get_module_repository_info,
    resolve_module_repository,
)
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

    with RepositoryManager(prefix="omnibenchmark_check") as repo_manager:
        for module_id, module in modules.items():
            logger.info(f"Validating module: {module_id}")
            modules_processed += 1

            # Get repository information
            repo_url, commit_hash = get_module_repository_info(benchmark_obj, module)

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

            # Resolve repository path (local or clone to temp)
            local_repo_path = resolve_module_repository(
                benchmark_obj, module, module_id, repo_manager
            )

            if local_repo_path is None:
                msg = f"Module {module_id}: Failed to access repository. Run the benchmark first to clone repositories."
                if warn:
                    validation_warnings.append(msg)
                    logger.warning(msg)
                else:
                    validation_errors.append(msg)
                    logger.error(msg)
                    ctx.exit(1)
                continue

            # Get file contents and presence information
            file_contents = repo_manager.get_repository_files(local_repo_path)
            files_present = repo_manager.get_files_present(local_repo_path)

            # Validate the module using existing validation functions
            try:
                result = validate_module_files(
                    module_id=module_id,
                    citation_content=file_contents.get("citation"),
                    license_content=file_contents.get("license"),
                    omnibenchmark_content=file_contents.get("omnibenchmark"),
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

    # Always cleanup at the end
    cleanup_temp_repositories("omnibenchmark_check")
