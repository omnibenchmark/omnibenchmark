"""Archive command for creating benchmark archives"""

import sys
import zipfile
from pathlib import Path

import click

from omnibenchmark.cli.remote import StorageAuth
from omnibenchmark.cli.utils.logging import logger
from omnibenchmark.remote.archive import archive_version
from omnibenchmark.remote.tree import tree_string_from_list

from .debug import add_debug_option


@add_debug_option
@click.command("archive")
@click.option(
    "--benchmark",
    "-b",
    required=True,
    type=click.Path(exists=True),
    help="Path to benchmark yaml file or benchmark id.",
    envvar="OB_BENCHMARK",
)
@click.option(
    "-c",
    "--code",
    help="Archive benchmarking code (repos).",
    is_flag=True,
    default=False,
    show_default=True,
)
@click.option(
    "-s",
    "--software",
    help="Archive software environments.",
    is_flag=True,
    default=False,
    show_default=True,
)
@click.option(
    "-r",
    "--results",
    help="Archive results files.",
    is_flag=True,
    default=False,
    show_default=True,
)
@click.option(
    "--compression",
    type=click.Choice(["none", "deflated", "bzip2", "lzma"], case_sensitive=False),
    default="none",
    help="Compression method.",
    show_default=True,
)
@click.option(
    "--compresslevel",
    type=int,
    default=None,
    help="Compression level.",
    show_default=True,
)
@click.option(
    "-n",
    "--dry-run",
    is_flag=True,
    default=False,
    help="Do not create the archive, just show what would be done.",
    show_default=True,
)
@click.option(
    "-r",
    "--use-remote-storage",
    help="Execute and store results remotely. Default False.",
    is_flag=True,
    default=False,
)
@click.option(
    "--out-dir", type=str, default="out", help="Output folder name (local only)."
)
@click.pass_context
def archive(
    ctx,
    benchmark,
    code,
    software,
    results,
    compression,
    compresslevel,
    dry_run,
    use_remote_storage,
    out_dir,
):
    """Archive a benchmark and its artifacts."""

    # Validate out_dir usage
    if use_remote_storage and out_dir != "out":
        logger.error("Invalid arguments: --out-dir can only be used with local storage")
        sys.exit(1)

    assert benchmark is not None

    storage_auth = StorageAuth(benchmark)
    benchmark = storage_auth.benchmark

    match compression:
        case "none":
            compression = zipfile.ZIP_STORED
        case "deflated":
            compression = zipfile.ZIP_DEFLATED
        case "bzip2":
            compression = zipfile.ZIP_BZIP2
        case "lzma":
            compression = zipfile.ZIP_LZMA
        case _:
            compression = zipfile.ZIP_STORED
    archive_file = archive_version(
        benchmark,
        outdir=Path("."),
        config=True,
        code=code,
        software=software,
        results=results,
        results_dir=out_dir,
        compression=compression,
        compresslevel=compresslevel,
        dry_run=dry_run,
        remote_storage=use_remote_storage,
    )
    if dry_run:
        click.echo(f"Files to archive:\n{tree_string_from_list(archive_file)}")
    else:
        click.echo(f"Created archive: {archive_file}")
