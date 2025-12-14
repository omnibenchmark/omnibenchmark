"""
Backward compatibility module for remote archive functionality.

This module maintains the original API while delegating to the new
archive module to avoid breaking existing code.
"""

# Import the new archive functionality
from omnibenchmark.archive import archive_benchmark

# Keep original imports for compatibility
import zipfile
from pathlib import Path
from typing import Optional


# Re-export the functions from the new module for backward compatibility
# These are already imported above


def archive_version(
    benchmark,
    outdir: Path = Path(),
    config: bool = True,
    code: bool = False,
    software: bool = False,
    results: bool = False,
    results_dir: str = "out",
    compression=zipfile.ZIP_STORED,
    compresslevel: Optional[int] = None,
    dry_run: bool = False,
    remote_storage: bool = True,
):
    """
    Legacy function name - delegates to archive_benchmark.

    Maintained for backward compatibility.
    """
    result = archive_benchmark(
        benchmark=benchmark,
        outdir=outdir,
        config=config,
        code=code,
        software=software,
        results=results,
        results_dir=results_dir,
        compression=compression,
        compresslevel=compresslevel,
        dry_run=dry_run,
        remote_storage=remote_storage,
    )

    if dry_run:
        return result
    else:
        return result[0] if result else None
