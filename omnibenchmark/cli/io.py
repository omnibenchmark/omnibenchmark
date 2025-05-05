"""cli commands related to input/output files"""

import json
import sys
import zipfile
from pathlib import Path

import click
import yaml

from omnibenchmark.benchmark import Benchmark
from omnibenchmark.cli.utils.logging import logger
from omnibenchmark.cli.utils.validation import validate_benchmark
from omnibenchmark.io.archive import archive_version
from omnibenchmark.io.files import checksum_files
from omnibenchmark.io.files import list_files
from omnibenchmark.io.files import download_files
from omnibenchmark.io.S3config import benchmarker_access_token_policy
from omnibenchmark.io.utils import tree_string_from_list
from omnibenchmark.io.utils import get_storage, remote_storage_args

from .debug import add_debug_option


@click.group(name="storage")
@click.pass_context
def storage(ctx):
    """Manage remote storage."""
    ctx.ensure_object(dict)


@add_debug_option
@storage.command(name="create-version")
@click.option(
    "-b",
    "--benchmark",
    help="Path to benchmark yaml file or benchmark id.",
    type=click.Path(exists=True),
    required=True,
    envvar="OB_BENCHMARK",
)
def create_benchmark_version(benchmark: str):
    """Create a new benchmark version."""

    with open(benchmark, "r") as fh:
        yaml.safe_load(fh)
        benchmark = Benchmark(Path(benchmark))

    auth_options = remote_storage_args(benchmark)
    # auth_options = {"endpoint": benchmark.converter.model.storage, "secure": False}

    # setup storage
    ss = get_storage(
        str(benchmark.converter.model.storage_api),
        auth_options,
        str(benchmark.converter.model.storage_bucket_name),
    )
    ss.set_version(benchmark.get_benchmark_version())
    if ss.version in ss.versions:
        logger.error(
            "Error: version already exists. Cannot overwrite.",
        )
        sys.exit(1)
    else:
        logger.info("Create a new benchmark version")
        ss.create_new_version(benchmark)


@add_debug_option
@storage.command(name="list")
@click.option(
    "-b",
    "--benchmark",
    help="Path to benchmark yaml file or benchmark id.",
    type=click.Path(exists=True),
    required=True,
    envvar="OB_BENCHMARK",
)
@click.option(
    "-t",
    "--type",
    help="File types. Options: all, code, inputs, outputs, logs, performance.",
    type=str,
    default="all",
)
@click.option("-s", "--stage", help="Stage to list files for.", type=str, default=None)
@click.option(
    "-i", "--id", "file_id", help="File id/type to list.", type=str, default=None
)
def list_all_files(
    benchmark: str,
    type: str = "all",
    stage: str = None,
    module: str = None,
    file_id: str = None,
):
    """List all or specific files for a benchmark."""
    if file_id is not None:
        logger.error("--file_id is not implemented")
        sys.exit(1)
    if type != "all":
        logger.error("--type is not implemented")
        sys.exit(1)

    objectnames, etags = list_files(
        benchmark=benchmark, type=type, stage=stage, module=module, file_id=file_id
    )
    if len(objectnames) > 0:
        for objectname, etag in zip(objectnames, etags):
            logger.info(f"{etag} {objectname}")


@add_debug_option
@storage.command(name="download")
@click.option(
    "-b",
    "--benchmark",
    help="Path to benchmark yaml file or benchmark id.",
    type=click.Path(exists=True),
    required=True,
    envvar="OB_BENCHMARK",
)
@click.option(
    "-t",
    "--type",
    help="File types. Options: all, code, inputs, outputs, logs, performance.",
    type=str,
    default="all",
)
@click.option(
    "-s", "--stage", help="Stage to download files from.", type=str, default=None
)
@click.option(
    "-m", "--module", help="Module to download files from.", type=str, default=None
)
@click.option(
    "-i", "--id", "file_id", help="File id to download.", type=str, default=None
)
@click.option(
    "-o",
    "--overwrite",
    help="Overwrite existing files.",
    is_flag=True,
    default=False,
    show_default=True,
)
def download_all_files(
    benchmark: str,
    type: str = "all",
    stage: str = None,
    module: str = None,
    file_id: str = None,
    overwrite: bool = False,
):
    """Download all or specific files for a benchmark."""
    if file_id is not None:
        logger.error("--file_id is not implemented")
        raise click.Abort()
    if type != "all":
        logger.info("--type is not implemented")
        raise click.Abort()

    download_files(
        benchmark=benchmark,
        type=type,
        stage=stage,
        module=module,
        file_id=file_id,
        verbose=True,
        overwrite=overwrite,
    )


@add_debug_option
@storage.command(name="checksum")
@click.option(
    "-b",
    "--benchmark",
    help="Path to benchmark yaml file or benchmark id.",
    type=click.Path(exists=True),
    required=True,
    envvar="OB_BENCHMARK",
)
def checksum_all_files(benchmark: str):
    """Generate md5sums of all benchmark outputs"""

    logger.info("Checking MD5 checksums... ")
    failed_checks_filenames = checksum_files(benchmark=benchmark, verbose=True)
    if len(failed_checks_filenames) > 0:
        logger.error("Failed checksums:")
        for filename in failed_checks_filenames:
            logger.error(filename)
        raise click.Abort()
    logger.info("Done")


@storage.command(name="create-policy")
@click.option(
    "-b",
    "--benchmark",
    help="Path to benchmark yaml file or benchmark id.",
    type=click.Path(exists=True),
    required=True,
    envvar="OB_BENCHMARK",
)
def create_policy(benchmark: str):
    """Create a new policy for a benchmark."""

    with open(benchmark, "r") as fh:
        yaml.safe_load(fh)
        benchmark = Benchmark(Path(benchmark))

    if (
        str(benchmark.converter.model.storage_api).upper() == "MINIO"
        or str(benchmark.converter.model.storage_api).upper() == "S3"
    ):
        policy = benchmarker_access_token_policy(
            benchmark.converter.model.storage_bucket_name
        )
        logger.error(json.dumps(policy, indent=2))
    else:
        logger.error("Error: Invalid storage type. Only MinIO/S3 storage is supported.")
        raise click.Abort()


@storage.command("archive")
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
    "-l",
    "--local",
    help="Execute and store results locally. Default False.",
    is_flag=True,
    default=False,
)
@click.pass_context
def archive_benchmark(
    ctx, benchmark, code, software, results, compression, compresslevel, dry_run, local
):
    """Archive a benchmark"""
    benchmark = validate_benchmark(benchmark, echo=False)

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
        compression=compression,
        compresslevel=compresslevel,
        dry_run=dry_run,
        local=local,
    )
    if dry_run:
        click.echo(f"Files to archive:\n{tree_string_from_list(archive_file)}")
    else:
        click.echo(f"Created archive: {archive_file}")
