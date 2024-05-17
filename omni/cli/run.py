"""cli commands related to benchmark/module execution and start"""

from pathlib import Path
from typing import List, Optional
from typing_extensions import Annotated

import typer

cli = typer.Typer(add_completion=False)


@cli.command("benchmark")
def run_benchmark(
    benchmark: Annotated[
        str,
        typer.Option(
            "--benchmark",
            "-b",
            help="Path to benchmark yaml file or benchmark id.",
        ),
    ],
    update: Annotated[
        bool,
        typer.Option(
            "--update",
            "-u",
            help="Run code for non existing outputs only.",
        ),
    ] = False,
    dry: Annotated[
        bool,
        typer.Option(
            "--dry",
            "-d",
            help="Dry run",
        ),
    ] = False,
    local: Annotated[
        bool,
        typer.Option(
            "--local",
            "-l",
            help="Execute and store results locally",
        ),
    ] = True,
    remote: Annotated[
        Optional[str],
        typer.Option(
            "--remote",
            "-r",
            help="The remote endpoint to use.",
        ),
    ] = None,
):
    """Run a benchmark as specified in the yaml."""
    typer.echo(f"Run {benchmark} in local {local}.", err=True)


@cli.command("module")
def run_module(
    benchmark: Annotated[
        str,
        typer.Option(
            "--benchmark",
            "-b",
            help="Path to benchmark yaml file or benchmark id.",
        ),
    ],
    repo: Annotated[
        str,
        typer.Option(
            "--repo",
            "-r",
            help="Repository url of the module to run",
        ),
    ],
    update: Annotated[
        bool,
        typer.Option(
            "--update",
            "-u",
            help="Run code for non existing outputs only.",
        ),
    ] = False,
    dry: Annotated[
        bool,
        typer.Option(
            "--dry",
            "-d",
            help="Dry run",
        ),
    ] = False,
    example: Annotated[
        bool,
        typer.Option(
            "--example",
            "-e",
            help="Run on example inputs only.",
        ),
    ] = True,
    all: Annotated[
        bool,
        typer.Option(
            "--all",
            "-a",
            help="Run on all valid benchmark inputs.",
        ),
    ] = False,
):
    """Run a specific module on all or example inputs locally."""
    typer.echo(
        f"Run {repo} as part of {benchmark} on example {example} inputs.", err=True
    )
    # NOTE: Do we also need a stage argument?
    # --all and --example are mutually exclusive. We use one flag only and run the other on default?
