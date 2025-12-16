"""CLI command for citation metadata extraction."""

from pathlib import Path

import click

from omnibenchmark.benchmark import BenchmarkExecution
from omnibenchmark.benchmark.cite import (
    extract_citation_metadata,
    format_output,
    CitationExtractionError,
)
from omnibenchmark.cli.utils.logging import logger


@click.command(name="cite")
@click.option(
    "--benchmark",
    "-b",
    required=True,
    type=click.Path(exists=True),
    help="Path to benchmark yaml file or benchmark id.",
    envvar="OB_BENCHMARK",
)
@click.option(
    "--format",
    "-f",
    type=click.Choice(["json", "yaml", "bibtex"], case_sensitive=False),
    default="yaml",
    help="Output format for citation information.",
)
@click.option(
    "--strict",
    is_flag=True,
    default=False,
    help="Fail on errors instead of converting them to warnings (default: warn mode).",
)
@click.option(
    "--out",
    type=click.Path(),
    help="Output file to write results to.",
)
@click.pass_context
def cite(ctx, benchmark: str, format: str, strict: bool, out: str):
    """Extract citation metadata from CITATION.cff files in benchmark modules.

    Automatically clones repositories to temporary directories if not found locally.
    By default, runs in warn mode and converts errors to warnings.
    Use --strict to fail on errors instead.
    """

    b = BenchmarkExecution(benchmark_yaml=Path(benchmark))
    if b is None:
        return

    logger.info(f"Extracting citation metadata from {benchmark}")
    logger.info("Will clone repositories to temporary directories if not found locally")

    try:
        citation_metadata = extract_citation_metadata(
            b, strict=strict, warn_mode=not strict
        )
    except RuntimeWarning as e:
        logger.warning(str(e))
        return
    except CitationExtractionError as e:
        logger.error(str(e))
        logger.error(f"Failed modules: {', '.join(e.failed_modules)}")
        # Log detailed error information
        for issue in e.issues:
            logger.error(f"  - {issue.message}")
        ctx.exit(1)

    try:
        output = format_output(citation_metadata, format, benchmark=b)
    except ValueError as e:
        logger.error(str(e))
        ctx.exit(1)

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
