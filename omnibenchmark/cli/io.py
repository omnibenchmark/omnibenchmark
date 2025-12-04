"""cli commands related to input/output files"""

import json
import sys
import zipfile

from pathlib import Path
from typing import TYPE_CHECKING

import click

from omnibenchmark.benchmark import BenchmarkExecution
from omnibenchmark.cli.utils.logging import logger
from omnibenchmark.io.archive import archive_version
from omnibenchmark.io.files import checksum_files
from omnibenchmark.io.files import list_files
from omnibenchmark.io.files import download_files
from omnibenchmark.io.tree import tree_string_from_list

if TYPE_CHECKING:
    from omnibenchmark.io.MinIOStorage import MinIOStorage

from .debug import add_debug_option


class StorageAuth:
    """Convenience class for handling storage authentication and validation."""

    def __init__(self, benchmark_path: str):
        from omnibenchmark.io.storage import remote_storage_args

        self.benchmark_path = benchmark_path
        self.benchmark = BenchmarkExecution(Path(benchmark_path))
        self.auth_options = remote_storage_args(self.benchmark)

        # Validate required storage components
        api = self.benchmark.get_storage_api()
        bucket = self.benchmark.get_storage_bucket_name()

        if api is None:
            logger.error("Error: No storage API found.")
            sys.exit(1)
        if bucket is None:
            logger.error("Error: No storage bucket found.")
            sys.exit(1)

        # Store validated non-null values
        self.api: str = api
        self.bucket: str = bucket

    def get_storage_instance(self) -> "MinIOStorage":
        """Get validated storage instance."""
        from omnibenchmark.io.storage import get_storage

        ss = get_storage(self.api, self.auth_options, self.bucket)
        if ss is None:
            logger.error("Error: No storage found.")
            sys.exit(1)
        # Type assertion since we know ss is not None after the exit check
        assert ss is not None
        return ss


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
    assert benchmark is not None

    storage_auth = StorageAuth(benchmark)
    ss = storage_auth.get_storage_instance()

    ss.set_version(storage_auth.benchmark.get_benchmark_version())

    if ss.version in ss.versions:
        logger.error(
            "Error: version already exists. Cannot overwrite.",
        )
        sys.exit(1)
    logger.info("Create a new benchmark version")
    ss.create_new_version(storage_auth.benchmark)


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
    stage: str = "",
    module: str = "",
    file_id: str = "",
):
    """List all or specific files for a benchmark."""
    if file_id is not None:
        logger.error("--file_id is not implemented")
        sys.exit(1)
    if type != "all":
        logger.error("--type is not implemented")
        sys.exit(1)

    objectnames, etags = list_files(benchmark, type, stage, module, file_id)
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
    stage: str = "",
    module: str = "",
    file_id: str = "",
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
        benchmark,
        type,
        stage,
        module,
        file_id,
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

    # ARCHITECTURAL NOTE: Business Logic in CLI Layer
    # This function contains domain logic (checksum validation) that belongs in
    # a service layer. During the LinkML → Pydantic migration, CLI commands retained
    # business logic to maintain functionality while models were refactored.
    #
    # Future consideration: Extract to a ChecksumService that can be used by CLI,
    # API endpoints, or other interfaces. CLI should only handle argument parsing
    # and output formatting.
    # TODO(ben): move this logic away from CLI
    logger.info("Checking MD5 checksums... ")
    failed_checks_filenames = checksum_files(
        benchmark=benchmark, type="all", stage="", module="", file_id="", verbose=True
    )
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
    "benchmark_path",
    help="Path to benchmark yaml file or benchmark id.",
    type=click.Path(exists=True),
    required=True,
    envvar="OB_BENCHMARK",
)
def create_policy(benchmark_path: str):
    """Create a new policy for a benchmark."""
    from omnibenchmark.io.S3config import benchmarker_access_token_policy

    assert benchmark_path is not None

    storage_auth = StorageAuth(benchmark_path)
    if storage_auth.api and (
        storage_auth.api.upper() == "MINIO" or storage_auth.api.upper() == "S3"
    ):
        policy = benchmarker_access_token_policy(storage_auth.bucket)
        logger.error(json.dumps(policy, indent=2))
    else:
        # TODO: this belongs to validation
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
    "--local-storage",
    help="Execute and store results locally. Default False.",
    is_flag=True,
    default=False,
)
@click.option(
    "--out-dir", type=str, default="out", help="Output folder name (local only)."
)
@click.pass_context
def archive_benchmark(
    ctx,
    benchmark,
    code,
    software,
    results,
    compression,
    compresslevel,
    dry_run,
    local_storage,
    out_dir,
):
    # Validate out_dir usage
    if not local_storage and out_dir != "out":
        logger.error(
            "-Invalid arguments: --out-dir can only be used with --local_storage"
        )
        sys.exit(1)

    """Archive a benchmark"""

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
        local_storage=local_storage,
    )
    if dry_run:
        click.echo(f"Files to archive:\n{tree_string_from_list(archive_file)}")
    else:
        click.echo(f"Created archive: {archive_file}")
