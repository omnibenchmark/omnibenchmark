"""cli commands related to input/output files"""

from typing import List, Optional
from typing_extensions import Annotated

import typer

cli = typer.Typer(
    add_completion=True,
    no_args_is_help=True,
    pretty_exceptions_short=False,
    rich_markup_mode=None,
    pretty_exceptions_enable=False,
    help="List, download and check input/output files.",
)


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
        f"List {type} files for {benchmark} at stage {stage} from module {module}:",
        err=True,
    )


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
        f"Download {type} files for {benchmark} at stage {stage} from module {module}",
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
