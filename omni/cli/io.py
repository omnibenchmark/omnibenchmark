"""cli commands related to input/output files"""

import json
from pathlib import Path
from typing import List, Optional

import typer
import yaml
from typing_extensions import Annotated

import omni.io.files
from omni.benchmark import Benchmark
from omni.io.utils import get_storage, remote_storage_args

cli = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    pretty_exceptions_short=False,
    rich_markup_mode=None,
    pretty_exceptions_enable=False,
    help="List, download and check input/output files.",
)


@cli.command("create-version")
def create_benchmark_version(
    benchmark: Annotated[
        str,
        typer.Option(
            "--benchmark",
            "-b",
            help="Path to benchmark yaml file or benchmark id.",
        ),
    ]
):
    """Create a new benchmark version."""
    with open(benchmark, "r") as fh:
        yaml.safe_load(fh)
        benchmark = Benchmark(Path(benchmark))

    auth_options = remote_storage_args(benchmark)
    # auth_options = {"endpoint": benchmark.converter.model.storage, "secure": False}

    # setup storage
    # TODO: add bucket_name to schema
    ss = get_storage(
        str(benchmark.converter.model.storage_api),
        auth_options,
        str(benchmark.converter.model.storage_bucket_name),
    )
    ss.set_version(benchmark.get_benchmark_version())
    if ss.version in ss.versions:
        typer.echo(
            "Error: version already exists. Cannot overwrite.",
            err=True,
        )
        raise typer.Exit(code=1)
    else:
        typer.echo("Create a new benchmark version")
        ss.create_new_version()


@cli.command("list")
def list_files(
    benchmark: Annotated[
        str,
        typer.Option(
            "--benchmark",
            "-b",
            help="Path to benchmark yaml file or benchmark id.",
        ),
    ],
    type: Annotated[
        str,
        typer.Option(
            "--type",
            "-t",
            help="File types. Options: all, code, inputs, outputs, logs, performance.",
        ),
    ] = "all",
    stage: Annotated[
        str,
        typer.Option(
            "--stage",
            "-s",
            help="Stage to list files for.",
        ),
    ] = None,
    module: Annotated[
        str,
        typer.Option(
            "--module",
            "-m",
            help="Module to list files for.",
        ),
    ] = None,
    file_id: Annotated[
        str,
        typer.Option(
            "--id",
            "-i",
            help="File id/type to list.",
        ),
    ] = None,
):
    """List all or specific files for a benchmark."""
    if file_id is not None:
        typer.echo("--file_id is not implemented")
        raise typer.Exit(code=1)
    if type != "all":
        typer.echo("--type is not implemented")
        raise typer.Exit(code=1)

    objectnames, etags = omni.io.files.list_files(
        benchmark=benchmark, type=type, stage=stage, module=module, file_id=file_id
    )
    if len(objectnames) > 0:
        for objectname, etag in zip(objectnames, etags):
            typer.echo(f"{etag} {objectname}")


@cli.command("download")
def download_files(
    benchmark: Annotated[
        str,
        typer.Option(
            "--benchmark",
            "-b",
            help="Path to benchmark yaml file or benchmark id.",
        ),
    ],
    type: Annotated[
        str,
        typer.Option(
            "--type",
            "-t",
            help="File types. Options: all, code, inputs, outputs, logs, performance.",
        ),
    ] = "all",
    stage: Annotated[
        str,
        typer.Option(
            "--stage",
            "-s",
            help="Stage to download files from.",
        ),
    ] = None,
    module: Annotated[
        Optional[str],
        typer.Option(
            "--module",
            "-m",
            help="Module to download files from.",
        ),
    ] = None,
    file_id: Annotated[
        Optional[List[str]],
        typer.Option(
            "--id",
            "-i",
            help="File id to download.",
        ),
    ] = None,
):
    """Download all or specific files for a benchmark."""
    if file_id is not None:
        typer.echo("--file_id is not implemented")
        raise typer.Exit(code=1)
    if type != "all":
        typer.echo("--type is not implemented")
        raise typer.Exit(code=1)

    omni.io.files.download_files(
        benchmark=benchmark,
        type=type,
        stage=stage,
        module=module,
        file_id=file_id,
        verbose=True,
    )


@cli.command("checksum")
def checksum_files(
    benchmark: Annotated[
        str,
        typer.Option(
            "--benchmark",
            "-b",
            help="Path to benchmark yaml file or benchmark id.",
        ),
    ],
):
    """Generate md5sums of all benchmark outputs"""
    typer.echo(f"Checking MD5 checksums... ", nl=False)
    failed_checks_filenames = omni.io.files.checksum_files(
        benchmark=benchmark, verbose=True
    )
    if len(failed_checks_filenames) > 0:
        typer.echo("Failed checksums:")
        for filename in failed_checks_filenames:
            typer.echo(filename)
        raise typer.Exit(code=1)
    typer.echo("Done")


@cli.command("create-policy")
def create_policy(
    benchmark: Annotated[
        str,
        typer.Option(
            "--benchmark",
            "-b",
            help="Path to benchmark yaml file or benchmark id.",
        ),
    ]
):
    """Create a new policy for a benchmark."""
    with open(benchmark, "r") as fh:
        yaml.safe_load(fh)
        benchmark = Benchmark(Path(benchmark))

    if (
        str(benchmark.converter.model.storage_api).upper() == "MINIO"
        or str(benchmark.converter.model.storage_api).upper() == "S3"
    ):
        policy = omni.io.S3config.benchmarker_access_token_policy(
            benchmark.converter.model.storage_bucket_name
        )
        typer.echo(json.dumps(policy, indent=2))
    else:
        typer.echo(
            "Error: Invalid storage type. Only MinIO/S3 storage is supported.", err=True
        )
        raise typer.Exit(code=1)
