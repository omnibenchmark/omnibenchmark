"""cli commands related to validation"""

from typing_extensions import Annotated

import typer

cli = typer.Typer(add_completion=False)


@cli.command("file")
def validate_file(
    benchmark: Annotated[
        str,
        typer.Option(
            "--benchmark",
            "-b",
            help="Path to benchmark yaml file or benchmark id.",
        ),
    ],
    stage: Annotated[
        str,
        typer.Option(
            "--stage",
            "-s",
            help="Stage to validate for",
        ),
    ],
    file_name: Annotated[
        str,
        typer.Option(
            "--filename",
            "-f",
            help="File name to check.",
        ),
    ],
    id: Annotated[
        str,
        typer.Option(
            "--id",
            "-i",
            help="File type id to check the file for.",
        ),
    ],
):
    """Validate file according to the benchmark stage and file type."""
    typer.echo(
        f"Validate {file_name} for {benchmark} stage {stage} and filetype {id}.",
        err=True,
    )


@cli.command("yaml")
def validate_yaml(
    benchmark: Annotated[
        str,
        typer.Option(
            "--benchmark",
            "-b",
            help="Path to benchmark yaml file or benchmark id.",
        ),
    ],
):
    """Validate benchmark yaml file."""
    typer.echo(f"Validate {benchmark}.", err=True)
