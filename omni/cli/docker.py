"""cli commands related to docker management"""

from typing_extensions import Annotated

import typer

cli = typer.Typer(add_completion=False)


@cli.command("download")
def download_docker(
    benchmark: Annotated[
        str,
        typer.Option(
            "--benchmark",
            "-b",
            help="Path to benchmark yaml file or benchmark id.",
        ),
    ],
):
    """Download a Docker image that fullfills the benchmarks software specifications.."""
    typer.echo(f"Download docker image for benchmark {benchmark}.")
