"""cli commands related to benchmark/module execution and start"""

from typing import Optional

import yaml
from typing_extensions import Annotated

import typer

from omni.benchmark import Benchmark
from omni.workflow.snakemake import SnakemakeEngine

cli = typer.Typer(add_completion=False)
workflow = SnakemakeEngine()


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
    if benchmark.endswith('.yaml') or benchmark.endswith('.yml'):
        try:
            with open(benchmark, 'r') as file:
                yaml.safe_load(file)
                benchmark = Benchmark(benchmark)
                typer.echo("Benchmark YAML file integrity check passed.")
                typer.echo("Running benchmark...")

                # TODO Include `local` and `remote` execution options once storage is integrated
                success = workflow.run_workflow(
                    benchmark,
                    cores=1, # TODO How should we configure this from the CLI? Should we leave the default to 1 core? What about other resources like memory?
                    update=update,
                    dry=dry
                )

                if success:
                    typer.echo("Benchmark run has finished successfully.", color=typer.colors.GREEN)
                else:
                    typer.echo("Benchmark run has failed.", err=True, color=typer.colors.RED)

                raise typer.Exit(code=success)

        except ValueError as e:
            typer.echo(f"Failed to parse YAML as a valid OmniBenchmark: {e}", err=True, color=typer.colors.RED)
            raise typer.Exit(code=1)

        except yaml.YAMLError as e:
            typer.echo(f"Error in YAML file: {e}.", err=True, color=typer.colors.RED)
            raise typer.Exit(code=1)

        except FileNotFoundError:
            typer.echo("Benchmark YAML file not found.", err=True, color=typer.colors.RED)
            raise typer.Exit(code=1)

        except Exception as e:
            typer.echo(f"An unexpected error occurred: {e}", err=True, color=typer.colors.RED)
            raise typer.Exit(code=1)

    else:
        typer.echo("Invalid benchmark input. Please provide a valid YAML file path.", err=True, color=typer.colors.RED)
        raise typer.Exit(code=1)


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
    stage: Annotated[
        str,
        typer.Option(
            "--stage",
            "-s",
            help="Stage to install software for.",
        ),
    ] = None,
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
    # --all and --example are mutually exclusive. Can we use one flag only and run the other on default?
