"""cli commands related to benchmark/module execution and start"""
import os
from pathlib import Path
from typing import Optional

import yaml
from typing_extensions import Annotated

import typer

from omni.benchmark import Benchmark
from omni.workflow.snakemake import SnakemakeEngine
from omni.workflow.workflow import WorkflowEngine

cli = typer.Typer(add_completion=False)
workflow: WorkflowEngine = SnakemakeEngine()


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
            help="The remote endpoint to use for storing results.",
        ),
    ] = None,
):
    """Run a benchmark as specified in the yaml."""

    # TODO Include `local` and `remote` execution options once storage is integrated
    # TODO Is it really necessary to have both `local` and `remote`. Could we assume that given an endpoint, you want the output files to be stored remotely?
    if not local:
        typer.echo(
            f"`local` argument is not supported yet. Workflows can only be run in local mode.",
            err=True,
        )
        raise typer.Exit(code=1)

    if remote is not None:
        typer.echo(
            f"`remote` argument is not supported yet. Workflows can only be run in local mode.",
            err=True,
        )
        raise typer.Exit(code=1)

    benchmark = validate_benchmark(benchmark)

    # TODO How should we configure `cores` from the CLI? Should we leave the default to 1 core? What about other resources like memory?
    typer.echo("Running benchmark...")
    success = workflow.run_workflow(benchmark, cores=1, update=update, dry=dry)

    if success:
        typer.echo("Benchmark run has finished successfully.", color=typer.colors.GREEN)
    else:
        typer.echo("Benchmark run has failed.", err=True, color=typer.colors.RED)

    raise typer.Exit(code=success)


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
    module: Annotated[
        str,
        typer.Option(
            "--module",
            "-m",
            help="Module (ID) to execute.",
        ),
    ],
    input_dir: Annotated[
        str,
        typer.Option(
            "--input",
            "-i",
            help="Input directory containing required input files for module.",
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

    # TODO Do we also need a stage argument?
    # TODO --all and --example are mutually exclusive. Can we use one flag only and run the other on default?
    if not example:
        typer.echo(
            f"`example` argument is not supported yet. Module workflows can only be run on example input.",
            err=True,
        )
        raise typer.Exit(code=1)

    if all is not None:
        typer.echo(
            f"`all` argument is not supported yet. Module workflows can only be run on example input.",
            err=True,
        )
        raise typer.Exit(code=1)

    benchmark = validate_benchmark(benchmark)
    benchmark_nodes = benchmark.get_nodes_by_module_id(module_id=module)
    if len(benchmark_nodes) > 0:
        typer.echo(
            f"Found {len(benchmark_nodes)} workflow nodes for module {module}. Running module benchmark..."
        )

        # Check available files in input_dir to figure out what dataset are we processing
        benchmark_datasets = benchmark.get_benchmark_datasets()
        dataset = None
        if os.path.isdir(input_dir):
            files = os.listdir(input_dir)
            base_names = [file.split(".")[0] for file in files]
            dataset = next((d for d in benchmark_datasets if d in base_names), None)

        if dataset is not None:
            for benchmark_node in benchmark_nodes:
                # TODO How should we configure `cores` from the CLI? Should we leave the default to 1 core? What about other resources like memory?
                success = workflow.run_node_workflow(
                    node=benchmark_node,
                    input_dir=input_dir,
                    dataset=dataset,
                    cores=1,
                    update=update,
                    dry=dry,
                )

                if success:
                    typer.echo(
                        "Module run has finished successfully.",
                        color=typer.colors.GREEN,
                    )
                else:
                    typer.echo(
                        "Module run has failed.", err=True, color=typer.colors.RED
                    )

                raise typer.Exit(code=success)
        else:
            typer.echo(
                f"Could not infer the appropriate dataset to run the node workflow on based on the files available in `{input_dir}`. None of the available datasets {benchmark_datasets} match the base names of the files.",
                err=True,
            )

            raise typer.Exit(code=1)

    else:
        typer.echo(
            f"Could not find module with id {module} in benchmark definition", err=True
        )
        raise typer.Exit(code=1)


def validate_benchmark(benchmark_file: str) -> Benchmark:
    if benchmark_file.endswith(".yaml") or benchmark_file.endswith(".yml"):
        try:
            with open(benchmark_file, "r") as file:
                yaml.safe_load(file)
                benchmark = Benchmark(Path(benchmark_file))
                typer.echo("Benchmark YAML file integrity check passed.")

                return benchmark

        except ValueError as e:
            typer.echo(
                f"Failed to parse YAML as a valid OmniBenchmark: {e}",
                err=True,
                color=typer.colors.RED,
            )
            raise typer.Exit(code=1)

        except yaml.YAMLError as e:
            typer.echo(f"Error in YAML file: {e}.", err=True, color=typer.colors.RED)
            raise typer.Exit(code=1)

        except FileNotFoundError:
            typer.echo(
                "Benchmark YAML file not found.", err=True, color=typer.colors.RED
            )
            raise typer.Exit(code=1)

        except Exception as e:
            typer.echo(
                f"An unexpected error occurred: {e}", err=True, color=typer.colors.RED
            )
            raise typer.Exit(code=1)

    else:
        typer.echo(
            "Invalid benchmark input. Please provide a valid YAML file path.",
            err=True,
            color=typer.colors.RED,
        )
        raise typer.Exit(code=1)
