"""cli commands related to input/output files"""

import json
from pathlib import Path

import click
import yaml


@click.group(name="storage")
@click.pass_context
def storage(ctx):
    """Manage remote storage."""
    ctx.ensure_object(dict)


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

    from omni.benchmark import Benchmark
    from omni.io.utils import get_storage, remote_storage_args

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
        click.echo(
            "Error: version already exists. Cannot overwrite.",
            err=True,
        )
        raise click.Exit(code=1)
    else:
        click.echo("Create a new benchmark version")
        ss.create_new_version()


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
@click.option("-i", "--id", help="File id/type to list.", type=str, default=None)
def list_files(
    benchmark: str,
    type: str = "all",
    stage: str = None,
    module: str = None,
    file_id: str = None,
):
    """List all or specific files for a benchmark."""
    if file_id is not None:
        click.echo("--file_id is not implemented")
        raise click.Exit(code=1)
    if type != "all":
        click.echo("--type is not implemented")
        raise click.Exit(code=1)

    objectnames, etags = omni.io.files.list_files(
        benchmark=benchmark, type=type, stage=stage, module=module, file_id=file_id
    )
    if len(objectnames) > 0:
        for objectname, etag in zip(objectnames, etags):
            click.echo(f"{etag} {objectname}")


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
@click.option("-i", "--id", help="File id to download.", type=str, default=None)
def download_files(
    benchmark: str,
    type: str = "all",
    stage: str = None,
    module: str = None,
    file_id: str = None,
):
    """Download all or specific files for a benchmark."""
    if file_id is not None:
        click.echo("--file_id is not implemented")
        raise click.Abort()
    if type != "all":
        click.echo("--type is not implemented")
        raise click.Abort()
    from omni.io.files import download_files

    download_files(
        benchmark=benchmark,
        type=type,
        stage=stage,
        module=module,
        file_id=file_id,
        verbose=True,
    )


@storage.command(name="checksum")
@click.option(
    "-b",
    "--benchmark",
    help="Path to benchmark yaml file or benchmark id.",
    type=click.Path(exists=True),
    required=True,
    envvar="OB_BENCHMARK",
)
def checksum_files(benchmark: str):
    """Generate md5sums of all benchmark outputs"""
    from omni.io.files import checksum_files

    click.echo(f"Checking MD5 checksums... ", nl=False)
    failed_checks_filenames = checksum_files(benchmark=benchmark, verbose=True)
    if len(failed_checks_filenames) > 0:
        click.echo("Failed checksums:")
        for filename in failed_checks_filenames:
            click.echo(filename)
        raise click.Abort()
    click.echo("Done")


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

    from omni.benchmark import Benchmark
    from omni.io.S3config import benchmarker_access_token_policy

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
        click.echo(json.dumps(policy, indent=2))
    else:
        click.echo(
            "Error: Invalid storage type. Only MinIO/S3 storage is supported.", err=True
        )
        raise click.Abort()
