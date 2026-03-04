"""CLI commands for generating benchmark dashboards"""

import csv
import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import click

from omnibenchmark.benchmark import BenchmarkExecution
from omnibenchmark.benchmark.dashboard import create_bettr_dashboard
from omnibenchmark.cli.utils.logging import logger
from omnibenchmark.cli.error_formatting import pretty_print_parse_error
from omnibenchmark.model.validation import BenchmarkParseError
from .debug import add_debug_option


DASHBOARD_FORMAT_EXT_DICT = {"bettr": "json"}


def _read_performance_tsv(file_path: Path) -> List[Dict[str, Any]]:
    """Read a TSV performance file as a list of dicts with numeric coercion."""
    rows: List[Dict[str, Any]] = []
    with open(file_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            parsed: Dict[str, Any] = {}
            for key, value in row.items():
                if value is None or value == "" or value.lower() == "nan":
                    parsed[key] = None
                else:
                    try:
                        parsed[key] = float(value)
                    except (ValueError, TypeError):
                        parsed[key] = value
            rows.append(parsed)
    return rows


@add_debug_option
@click.command(name="dashboard")
@click.argument("benchmark", type=click.Path(exists=True))
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

    logger.info(
        f"Generating {dashboard_format} dashboard from benchmark performance results..."
    )

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

    try:
        dashboard_data = {}
        rows = _read_performance_tsv(performance_file)
        if dashboard_format.lower() == "bettr":
            dashboard_data = create_bettr_dashboard(rows, id_col="module")
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
