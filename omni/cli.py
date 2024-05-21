"""Cli implementation of omni-py via typer"""

from pathlib import Path
from typing import Any, List, Mapping, Optional
from typing_extensions import Annotated

import typer

cli = typer.Typer(add_completion=False)


@cli.command()
def fetch(
    module_url: Annotated[str, typer.Argument(...)],
    commit: Annotated[
        str,
        typer.Option(
            "--commit",
            "-v",
            help="Remote commit to fetch from. Defaults to the latest commit on the default branch.",
        ),
    ] = "latest",
):
    """Fetch a specific version of a benchmark module"""
    typer.echo(f"Fetch module with url {module_url}", err=True)
    typer.echo(f"Using commit: {commit}", err=True)


@cli.command()
def list_inputs(
    stage_id: Annotated[str, typer.Argument(...)],
    benchmark: Annotated[
        str,
        typer.Argument(
            ...,
            # Question: Should this be the benchmark url?
            help="The benchmark name the satge belongs to",
        ),
    ],
    version: Annotated[
        str,
        typer.Option(
            "--version",
            "-v",
            help="Benchmark version to use",
        ),
    ] = "latest",
    test: Annotated[
        Optional[bool],
        typer.Option(
            "--test",
            "-t",
            help="Show test inputs only",
        ),
    ] = True,
    download: Annotated[
        str,
        typer.Option(
            "--download",
            "-d",
            help="Download files from remote.",
        ),
    ] = False,
):
    """List and/or download (test) input files for a specific stage and benchmark"""
    typer.echo(
        f"Get inputs for stage {stage_id} in benchmark {benchmark} {version}", err=True
    )


@cli.command()
def list_outputs(
    module_id: Annotated[str, typer.Argument(...)],
    benchmark: Annotated[
        str,
        typer.Argument(
            ...,
            # Question: Should this be the benchmark url?
            help="The benchmark name the satge belongs to",
        ),
    ],
    version: Annotated[
        str,
        typer.Option(
            "--version",
            "-v",
            help="Benchmark version to use",
        ),
    ] = "latest",
    download: Annotated[
        str,
        typer.Option(
            "--download",
            "-d",
            help="Download files from remote.",
        ),
    ] = False,
):
    """List or download (test) input files for a specific stage and benchmark"""
    typer.echo(
        f"Get outputs for module {module_id} in benchmark {benchmark} {version}",
        err=True,
    )


@cli.command()
def get(
    remote_paths: Annotated[List[Path], typer.Argument(...)],
):
    """Download explicit files from remote."""
    typer.echo(f"Get files {[fi for fi in remote_paths]}", err=True)


@cli.command()
def run(
    # Question: Should this be a path to a script instead?
    module_path: Annotated[Path, typer.Argument(...)],
    outputs: Annotated[
        List[str],
        typer.Option(
            "--output",
            "-o",
            help="The matched output type and paths, e.g. output1=path/to/output1.txt",
        ),
    ],
    inputs: Annotated[
        Optional[List[str]],
        typer.Option(
            "--input",
            "-i",
            help="The matched input type and paths, e.g. input1=path/to/input1.txt",
        ),
    ] = None,
    params: Annotated[
        Optional[List[str]],
        typer.Option(
            "--param",
            "-p",
            help="The matched output type and paths, e.g. param1=X",
        ),
    ] = None,
):
    """Execute a modules code with explicit inputs, outputs and parameter"""
    # TODO: Include the environment here?
    typer.echo(f"Run module {module_path}.", err=True)
    typer.echo(f"Outputs: {[out for out in outputs]}.", err=True)
    if inputs:
        typer.echo(f"Inputs: {[input for input in inputs]}.", err=True)
    if params:
        typer.echo(f"Params: {[param for param in params]}.", err=True)


@cli.command()
def upload(
    files: Annotated[list[Path], typer.Argument(...)],
    s3_endpoint: Annotated[
        Optional[str],
        typer.Option(
            "--s3_endpoint",
            "-s",
            help="S3 endpoint including bucket name to upload files to.",
        ),
    ] = None,
):
    """Upload one or more files to a S3 bucket."""
    typer.echo(f"Upload {[fi for fi in files]} to {s3_endpoint}", err=True)
