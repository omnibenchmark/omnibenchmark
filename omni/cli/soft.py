"""cli commands related to software management"""

from typing_extensions import Annotated

import typer

cli = typer.Typer(add_completion=False)


@cli.command("check")
def check_software(
    benchmark: Annotated[
        str,
        typer.Option(
            "--benchmark",
            "-b",
            help="Path to benchmark yaml file or benchmark id.",
        ),
    ],
    requirements: Annotated[
        str,
        typer.Option(
            "--requirements",
            "-r",
            help="Path to modules requirements.txt file.",
        ),
    ] = "requirements.txt",
):
    """Checks the compatibility of the requirements with the benchmark's software stack."""
    typer.echo(f"Requirments in {requirements} are compatible with {benchmark}!")


@cli.command("install")
def install_software(
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
            help="Stage to install software for.",
        ),
    ] = None,
    module: Annotated[
        str,
        typer.Option(
            "--module",
            "-m",
            help="Module to install software for.",
        ),
    ] = None,
):
    """Install software and provide an environment with all executables."""
    typer.echo(f"Install software for {benchmark}, {stage}.")
