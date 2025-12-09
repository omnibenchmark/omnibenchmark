"""CLI commands for generating benchmark dashboards"""

import json
import sys
from pathlib import Path

import click

from omnibenchmark.benchmark import BenchmarkExecution
from omnibenchmark.benchmark.dashboard import create_bettr_dashboard, EXPORT_AVAILABLE
from omnibenchmark.cli.utils.logging import logger
from omnibenchmark.cli.error_formatting import pretty_print_parse_error
from omnibenchmark.model.validation import BenchmarkParseError
from .debug import add_debug_option


DASHBOARD_FORMAT_EXT_DICT = {"bettr": "json"}


@add_debug_option
@click.command(name="dashboard")
@click.option(
    "-b",
    "--benchmark",
    help="Path to benchmark yaml file or benchmark id.",
    required=True,
    envvar="OB_BENCHMARK",
    type=click.Path(exists=True),
)
@click.option(
    "--format",
    "-f",
    "dashboard_format",
    type=click.Choice(["bettr"], case_sensitive=False),
    default="bettr",
    help="Dashboard format to generate (default: bettr).",
)
@click.option(
    "--out-dir",
    "-o",
    type=click.Path(),
    default="out",
    help="Output directory containing benchmark results (default: out).",
)
@click.pass_context
def dashboard(ctx, benchmark, dashboard_format, out_dir):
    """Generate a dashboard from benchmark results.

    This command generates a dashboard from the performance data
    collected during benchmark execution.
    """
    ctx.ensure_object(dict)

    if not EXPORT_AVAILABLE:
        logger.error(
            "Export functionality not available. Install with: pip install omnibenchmark[export]"
        )
        sys.exit(1)

    logger.info(
        f"Generating {dashboard_format} dashboard from benchmark performance results..."
    )

    import pandas

    try:
        _ = BenchmarkExecution(Path(benchmark), Path(out_dir))
        logger.info("Benchmark YAML file integrity check passed.")
    except BenchmarkParseError as e:
        formatted_error = pretty_print_parse_error(e)
        logger.error(f"Failed to load benchmark: {formatted_error}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Failed to load benchmark: {e}")
        sys.exit(1)

    out_path = Path(out_dir)

    if not out_path.exists():
        logger.error(f"Output directory does not exist: {out_dir}")
        sys.exit(1)

    performance_file = out_path / "performances.tsv"
    if not performance_file.exists():
        logger.error(
            f"No performance data found at {performance_file}. "
            "Please run the benchmark first to generate performance metrics."
        )
        sys.exit(1)

    # Generate dashboard based on format
    try:
        dashboard_data = {}
        performance_df = pandas.read_csv(performance_file, sep="\t")
        if dashboard_format.lower() == "bettr":
            dashboard_data = create_bettr_dashboard(performance_df, id_col="module")
        else:
            logger.error(f"Unsupported dashboard format: {dashboard_format}")
            sys.exit(1)

        output_ext = DASHBOARD_FORMAT_EXT_DICT[dashboard_format]
        output = out_path / f"{dashboard_format}_dashboard.{output_ext}"
        with open(output, "w") as f:
            json.dump(dashboard_data, f, indent=2)

        logger.info(f"Successfully generated {dashboard_format} dashboard.")
        logger.info(f"Dashboard saved to: {output}")
        sys.exit(0)

    except Exception as e:
        logger.error(f"Dashboard generation failed: {e}")
        sys.exit(1)
