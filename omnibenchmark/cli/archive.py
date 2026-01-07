"""Archive command for creating benchmark archives"""

import click
import zipfile
from pathlib import Path

from omnibenchmark.benchmark import BenchmarkExecution
from omnibenchmark.cli.remote import StorageAuth
from omnibenchmark.archive import archive_benchmark
from omnibenchmark.remote.tree import tree_string_from_list

from .debug import add_debug_option


@add_debug_option
@click.command("archive")
@click.argument(
    "benchmark",
    type=click.Path(exists=True),
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
    default="bzip2",
    help="Compression method.",
    show_default=True,
)
@click.option(
    "--compresslevel",
    type=int,
    default=9,
    help="Compression level (1-9, higher is better compression).",
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
    "--use-remote-storage",
    help="Download results from remote storage for archiving. Default False.",
    is_flag=True,
    default=False,
)
@click.option(
    "--out-dir",
    help="Local directory name to archive from (local-only mode). Default: `out`",
    default=None,
    type=str,
)
@click.option(
    "-o",
    "--output-file",
    help="Custom output file path for the archive. Extension must match compression type: none/deflated (.zip), bzip2 (.bz2), lzma (.xz).",
    default=None,
    type=click.Path(),
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
    output_file,
):
    """Archive a benchmark and its artifacts.

    BENCHMARK: Path to benchmark YAML file.
    """

    # Validate compression=none with compresslevel
    if compression == "none" and compresslevel != 9:
        raise click.UsageError(
            "Cannot specify --compresslevel with --compression=none. "
            "Remove --compresslevel or choose a different compression method."
        )

    # Validate out_dir usage
    if use_remote_storage and out_dir:
        raise click.UsageError("--out-dir can only be used with local storage")

    # Validate output_file extension against compression
    if output_file:
        output_path = Path(output_file)
        expected_extensions = {
            zipfile.ZIP_STORED: [".zip"],
            zipfile.ZIP_DEFLATED: [".zip"],
            zipfile.ZIP_BZIP2: [".bz2"],
            zipfile.ZIP_LZMA: [".xz"],
        }

        # Convert compression string to zipfile constant for validation
        compression_map = {
            "none": zipfile.ZIP_STORED,
            "deflated": zipfile.ZIP_DEFLATED,
            "bzip2": zipfile.ZIP_BZIP2,
            "lzma": zipfile.ZIP_LZMA,
        }

        compression_type = compression_map.get(compression, zipfile.ZIP_STORED)
        valid_extensions = expected_extensions.get(compression_type, [".zip"])

        if output_path.suffix.lower() not in valid_extensions:
            all_combinations = ["none/deflated → .zip", "bzip2 → .bz2", "lzma → .xz"]
            raise click.UsageError(
                f"Output file extension '{output_path.suffix}' does not match "
                f"compression type '{compression}'. Expected for '{compression}': {', '.join(valid_extensions)}.\n"
                f"Valid combinations: {' | '.join(all_combinations)}"
            )

    # Load benchmark differently based on whether we need storage authentication
    if use_remote_storage:
        # For remote storage, we need full storage authentication
        storage_auth = StorageAuth(benchmark)
        benchmark_obj = storage_auth.benchmark
    else:
        # For local-only archiving, we can load benchmark without storage requirements
        benchmark_obj = BenchmarkExecution(Path(benchmark))

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

    results_dir = out_dir if out_dir else "out"

    # Determine output directory and filename
    if output_file:
        output_path = Path(output_file)
        outdir = output_path.parent if output_path.parent != Path() else Path(".")
        custom_filename = output_path.name
    else:
        outdir = Path(".")
        custom_filename = None

    archive_file = archive_benchmark(
        benchmark_obj,
        outdir=outdir,
        config=True,
        code=code,
        software=software,
        results=results,
        results_dir=results_dir,
        compression=compression,
        compresslevel=compresslevel,
        dry_run=dry_run,
        remote_storage=use_remote_storage,
        custom_filename=custom_filename,
    )
    if dry_run:
        click.echo(f"Files to archive:\n{tree_string_from_list(archive_file)}")
    else:
        click.echo(f"Created archive: {archive_file[0]}")
