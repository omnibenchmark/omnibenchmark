"""
Script to aggregate performance files from all benchmark runs.
This is called by the Snakemake 'all' rule when aggregate_performance=True.

NOTE: This uses the `script:` directive instead of `run:` block to ensure
the code executes in the correct Python environment (the one running Snakemake).
"""

from pathlib import Path
import os
from typing import TYPE_CHECKING

from omnibenchmark.workflow.snakemake.scripts.parse_performance import (
    write_combined_performance_file,
)

if TYPE_CHECKING:
    pass


def main(snakemake):  # type: ignore[no-untyped-def]
    """Aggregate performance files.

    Args:
        snakemake: Snakemake object with input, output, wildcards, etc.
    """
    # Get the output directory from the first output file
    output_files = snakemake.output
    output_dir = Path(str(os.path.commonpath(output_files)))
    if len(output_files) == 1:
        output_dir = Path(os.path.dirname(output_dir))

    # Get the benchmark output directory from config or infer from paths
    out_dir = snakemake.config.get("out_dir", "out")

    # Find all performance files
    from snakemake.io import glob_wildcards, expand

    result = glob_wildcards(f"{out_dir}/{{path}}/{{dataset}}_performance.txt")
    performances_raw = expand(
        f"{out_dir}/{{path}}/{{dataset}}_performance.txt",
        path=result.path,
        dataset=result.dataset,
    )
    performances = sorted(list(set(performances_raw)))  # type: ignore[arg-type]

    # Write combined performance file
    write_combined_performance_file(output_dir, performances)


if __name__ == "__main__":
    # When used as a Snakemake script, snakemake object is automatically available
    main(snakemake)  # type: ignore[name-defined]  # noqa: F821
