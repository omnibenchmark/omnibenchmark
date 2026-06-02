"""Metric-collector output path formatting."""

from pathlib import Path


def format_mc_output(output, out_dir: Path, collector_id: str):
    """Format metric collector output path.

    Args:
        output: IOFile object
        out_dir: Output directory path
        collector_id: Collector identifier
    """
    if output.path:
        o = output.path.replace("{name}", collector_id)
        return f"{out_dir}/{o}"
    else:
        return str(out_dir / output.id)
