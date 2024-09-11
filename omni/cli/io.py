"""cli commands related to input/output files"""

from typing import List, Optional

import typer
from typing_extensions import Annotated

from omni.cli.run import validate_benchmark
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
    benchmark = validate_benchmark(benchmark)

    auth_options = remote_storage_args(benchmark)
    # auth_options = {"endpoint": benchmark.converter.model.storage, "secure": False}

    # setup storage
    # TODO: add bucket_name to schema
    ss = get_storage(
        str(benchmark.converter.model.storage_api),
        auth_options,
        str(benchmark.converter.model.id),
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
    typer.echo(
        f"(not implemented) List {type} files for {benchmark} at stage {stage} from module {module}:",
        err=True,
    )
    if file_id is not None:
        typer.echo("--file_id is not implemented")
        raise typer.Exit(code=1)
    if type is not None:
        typer.echo("--type is not implemented")
        raise typer.Exit(code=1)

    benchmark = validate_benchmark(benchmark)

    auth_options = {"endpoint": benchmark.converter.model.storage, "secure": False}

    # TODO: add bucket_name to schema
    ss = get_storage(
        str(benchmark.converter.model.storage_api),
        auth_options,
        str(benchmark.converter.model.id),
    )
    ss.set_version(benchmark.get_benchmark_version())
    ss._get_objects()

    all_files = benchmark.get_output_paths()
    expected_files = []
    for file in all_files:
        filter_stage = False
        if stage is not None:
            if file.split("/")[-4] != stage:
                filter_stage = True

        filter_module = False
        if module is not None:
            if file.split("/")[-3] != module:
                filter_module = True

        if not filter_stage and not filter_module:
            expected_files.append(file.replace("{dataset}", file.split("/")[2]))

    real_files = {k: v for k, v in ss.files.items() if k in expected_files}

    if len(real_files) > 0:
        for file in real_files:
            typer.echo(file)


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
    typer.echo(
        f" (not implemented) Download {type} files for {benchmark} at stage {stage} from module {module}",
        err=True,
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
    typer.echo(f"Generate md5sums for benchmark {benchmark}")
