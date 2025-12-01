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
    "--warn",
    is_flag=True,
    help="Convert errors to warnings instead of failing fast (default: strict mode fails on first error).",
)
@click.option(
    "--out",
    type=click.Path(),
    help="Output file to write results to.",
)
@click.pass_context
def cite(ctx, benchmark: str, format: str, warn: bool, out: str):
    """Extract citation metadata from CITATION.cff files in benchmark modules.

    Automatically clones repositories to temporary directories if not found locally.
    By default, runs in strict mode and fails fast on the first error found.
    Use --warn to convert errors to warnings and continue processing all modules.
    """

    b = BenchmarkExecution(benchmark_yaml=Path(benchmark))
    if b is None:
        return

    logger.info(f"Extracting citation metadata from {benchmark}")
    logger.info("Will clone repositories to temporary directories if not found locally")

    try:
        citation_metadata = extract_citation_metadata(b, strict=True, warn_mode=warn)
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
        output = format_output(citation_metadata, format)
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
