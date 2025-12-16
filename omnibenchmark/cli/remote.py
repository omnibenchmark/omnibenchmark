"""cli commands related to input/output files"""

import json
import sys

from pathlib import Path

import click

from omnibenchmark.benchmark import BenchmarkExecution
from omnibenchmark.cli.utils.logging import logger
from omnibenchmark.cli.error_formatting import pretty_print_parse_error
from omnibenchmark.model.validation import BenchmarkParseError
from omnibenchmark.remote.files import checksum_files
from omnibenchmark.remote.files import list_files
from omnibenchmark.remote.files import download_files
from datetime import datetime
from difflib import unified_diff

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from omnibenchmark.remote.MinIOStorage import MinIOStorage

from .debug import add_debug_option


class StorageAuth:
    """Convenience class for handling storage authentication and validation."""

    def __init__(self, benchmark_path: str, require_credentials: bool = True):
        from omnibenchmark.remote.storage import remote_storage_args

        self.benchmark_path = benchmark_path
        try:
            self.benchmark = BenchmarkExecution(Path(benchmark_path))
        except BenchmarkParseError as e:
            formatted_error = pretty_print_parse_error(e)
            logger.error(f"Failed to load benchmark:\n{formatted_error}")
            sys.exit(1)
        self.auth_options = remote_storage_args(
            self.benchmark, required=require_credentials
        )

        # Validate required storage components
        import click

        api = self.benchmark.get_storage_api()
        bucket = self.benchmark.get_storage_bucket_name()

        if api is None:
            logger.error(
                click.style("[ERROR]", fg="red", bold=True)
                + " No storage API configured. Set 'storage.api' in your benchmark YAML."
            )
            sys.exit(1)
        if bucket is None:
            logger.error(
                click.style("[ERROR]", fg="red", bold=True)
                + " No storage bucket configured. Set 'storage.bucket_name' in your benchmark YAML."
            )
            sys.exit(1)

        # Store validated non-null values
        self.api: str = api
        self.bucket: str = bucket

    def get_storage_instance(self) -> "MinIOStorage":
        """Get validated storage instance."""
        from omnibenchmark.remote.storage import get_storage

        ss = get_storage(self.api, self.auth_options, self.bucket)
        if ss is None:
            logger.error("Error: No storage found.")
            sys.exit(1)
        # Type assertion since we know ss is not None after the exit check
        assert ss is not None
        return ss


@click.group(name="remote")
@click.pass_context
def remote(ctx):
    """Manage remote storage."""
    ctx.ensure_object(dict)


@click.group(name="files")
@click.pass_context
def files(ctx):
    """Manage files in remote storage."""
    ctx.ensure_object(dict)


remote.add_command(files)


@click.group(name="version")
@click.pass_context
def version(ctx):
    """Manage benchmark versions."""
    ctx.ensure_object(dict)


remote.add_command(version)


@click.group(name="policy")
@click.pass_context
def policy(ctx):
    """Manage storage policies."""
    ctx.ensure_object(dict)


remote.add_command(policy)


@add_debug_option
@version.command(name="create")
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
@files.command(name="list")
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
@files.command(name="download")
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
@files.command(name="checksum")
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
    # TODO: Add CLI option to specify custom local directory (currently hardcoded to "out")
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


@add_debug_option
@policy.command(name="create")
@click.option(
    "-b",
    "--benchmark",
    "benchmark_path",
    help="Path to benchmark yaml file.",
    type=click.Path(exists=True),
    required=False,
    envvar="OB_BENCHMARK",
)
@click.option(
    "--bucket",
    "bucket_name",
    help="S3 bucket name. If not provided, reads from benchmark YAML.",
    type=str,
    required=False,
)
def create_policy(benchmark_path: str, bucket_name: str):
    """Generate an S3/MinIO IAM policy for a benchmark bucket.

    This command generates a least-privilege AWS IAM policy JSON that can be used
    to create access keys in MinIO or AWS. The policy allows full S3 operations
    on the bucket but denies deletion and governance bypass.

    You can either:
    - Provide a benchmark YAML file (-b) to read the bucket name from storage.bucket_name
    - Provide the bucket name directly (--bucket)
    """
    from omnibenchmark.remote.S3config import benchmarker_access_token_policy
    import click

    # Get bucket name from either parameter or benchmark YAML
    if bucket_name:
        # Direct bucket name provided
        bucket = bucket_name
    elif benchmark_path:
        # Load from benchmark YAML
        from omnibenchmark.model.benchmark import Benchmark

        try:
            benchmark = Benchmark.from_yaml(Path(benchmark_path))
        except BenchmarkParseError as e:
            formatted_error = pretty_print_parse_error(e)
            logger.error(f"Failed to load benchmark:\n{formatted_error}")
            sys.exit(1)

        # Get storage configuration
        api = benchmark.get_storage_api()
        bucket = benchmark.get_storage_bucket_name()

        # Validate we have bucket name
        if bucket is None:
            logger.error(
                click.style("[ERROR]", fg="red", bold=True)
                + " No storage bucket configured. Set 'storage.bucket_name' in your benchmark YAML or use --bucket."
            )
            sys.exit(1)

        # Validate storage API is S3/MinIO (both use AWS IAM policies)
        if api is None:
            logger.error(
                click.style("[ERROR]", fg="red", bold=True)
                + " Storage API not configured. Set 'storage.api' to 'S3' or 'MinIO' in your benchmark YAML."
            )
            sys.exit(1)

        if api.upper() not in ("MINIO", "S3"):
            logger.error(
                click.style("[ERROR]", fg="red", bold=True)
                + f" Storage API '{api}' does not use AWS IAM policies. Only S3 and MinIO are supported."
            )
            sys.exit(1)
    else:
        # Neither provided
        logger.error(
            click.style("[ERROR]", fg="red", bold=True)
            + " Must provide either --benchmark or --bucket option."
        )
        sys.exit(1)

    # Generate and print the policy
    policy = benchmarker_access_token_policy(bucket)
    print(json.dumps(policy, indent=2))


@add_debug_option
@version.command("diff")
@click.option(
    "--benchmark",
    "-b",
    "benchmark_path",
    required=True,
    type=click.Path(exists=True),
    help="Path to benchmark yaml file or benchmark id.",
    envvar="OB_BENCHMARK",
)
@click.option(
    "--version1",
    "-v1",
    required=True,
    type=str,
    help="Reference version.",
)
@click.option(
    "--version2",
    "-v2",
    required=True,
    type=str,
    help="Version to compare with.",
)
@click.pass_context
def diff_benchmark(ctx, benchmark_path: str, version1, version2):
    """Show differences between 2 benchmark versions."""
    from omnibenchmark.remote.storage import get_storage, remote_storage_args

    logger.info(
        f"Found the following differences in {benchmark_path} for {version1} and {version2}."
    )
    b = BenchmarkExecution(Path(benchmark_path))
    auth_options = remote_storage_args(benchmark_path)

    api = b.get_storage_api()
    bucket = b.get_storage_bucket_name()
    if api is None or bucket is None:
        raise (ValueError)
    # setup storage
    #
    # TODO: use walrus
    ss = get_storage(api, auth_options, bucket)
    if ss is not None:
        # get objects for first version
        ss.set_version(version1)
        ss._get_objects()
        files_v1 = [
            f"{f[0]}   {f[1]['size']}   {datetime.fromisoformat(f[1]['last_modified']).strftime('%Y-%m-%d %H:%M:%S')}\n"
            for f in ss.files.items()
        ]
        if f"versions/{version1}.csv" in ss.files.keys():
            creation_time_v1 = datetime.fromisoformat(
                ss.files[f"versions/{version1}.csv"]["last_modified"]
            ).strftime("%Y-%m-%d %H:%M:%S")
        else:
            creation_time_v1 = ""

        # get objects for second version
        ss.set_version(version2)
        ss._get_objects()
        files_v2 = [
            f"{f[0]}   {f[1]['size']}   {datetime.fromisoformat(f[1]['last_modified']).strftime('%Y-%m-%d %H:%M:%S')}\n"
            for f in ss.files.items()
        ]
        if f"versions/{version2}.csv" in ss.files.keys():
            creation_time_v2 = datetime.fromisoformat(
                ss.files[f"versions/{version2}.csv"]["last_modified"]
            ).strftime("%Y-%m-%d %H:%M:%S")
        else:
            creation_time_v2 = ""

    # diff the two versions
    click.echo(
        "".join(
            list(
                unified_diff(
                    files_v1,
                    files_v2,
                    fromfile=f"version {version1}",
                    tofile=f"version {version2}",
                    fromfiledate=creation_time_v1,
                    tofiledate=creation_time_v2,
                )
            )
        )
    )


@add_debug_option
@version.command("list")
@click.option(
    "--benchmark",
    "-b",
    "benchmark_path",
    required=True,
    type=click.Path(exists=True),
    help="Path to benchmark yaml file or benchmark id.",
    envvar="OB_BENCHMARK",
)
@click.pass_context
def list_versions(ctx, benchmark_path: str):
    """List all available benchmark versions."""
    from omnibenchmark.remote.storage import get_storage, remote_storage_args

    logger.info(f"Available versions of {benchmark_path}:")

    b = BenchmarkExecution(Path(benchmark_path))
    auth_options = remote_storage_args(b)
    api = b.get_storage_api()
    bucket = b.get_storage_bucket_name()

    if api is None:
        raise ValueError("No storage API found")
    if bucket is None:
        raise ValueError("No storage bucket found")

    # setup storage
    ss = get_storage(api, auth_options, bucket)
    if ss is None:
        raise ValueError("No storage found")

    if len(ss.versions) > 0:
        if len(ss.versions) > 1:
            ss.versions.sort()
        for version in ss.versions:
            click.echo(f"{str(version):>8}")
