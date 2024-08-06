"""Cli implementation of omni-py via typer"""

import typer

from pathlib import Path
from typing_extensions import Annotated

from omni.cli import benchmark
from omni.cli import docker
from omni.cli import soft
from omni.cli import io
from omni.cli import run
from omni.cli import validate


cli = typer.Typer(add_completion=False, no_args_is_help=True)
cli.add_typer(
    benchmark.cli, name="benchmark", help="Manage benchmarks and their versions."
)
cli.add_typer(
    soft.cli,
    name="software",
    help="Manage and install benchmark-specific software environments",
)
# cli.add_typer(
#     docker.cli, name="docker", help="use docker to manage software environments"
# )
cli.add_typer(io.cli, name="files", help="List, download and check input/output files.")
cli.add_typer(run.cli, name="run", help="Execute benchmarks or modules")
cli.add_typer(
    validate.cli, name="validate", help="Validate benchmarks, modules or files"
)


@cli.command("trigger checks")
def trigger_checks(
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
            help="Repository url to check.",
        ),
    ],
):
    """Trigger benchmark-specific checks in a certain repository."""
    typer.echo(f"Trigger {benchmark} checks in {repo}.", err=True)


@cli.command("start module")
def start_module(
    module_name: Annotated[Path, typer.Argument(...)],
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
            help="Stage to associate the module with.",
        ),
    ],
):
    """Start a new gitlab module at a certain stage."""
    typer.echo(
        f"Start module {module_name} as part of {benchmark} on stage {stage}.", err=True
    )
    # NOTE: We probably also need a gitlab url? Not sure about module_name?
