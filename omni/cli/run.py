"""cli commands related to benchmark/module execution and start"""
import os
from itertools import chain
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
    cores: Annotated[
        int,
        typer.Option(
            "--cores",
            "-c",
            help="The parallelism level for the workflow scheduler.",
        ),
    ] = 1,
    update: Annotated[
        bool,
        typer.Option(
            "--update",
            "-u",
            help="Force re-run execution for all modules and stages.",
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
            help="Execute and store results locally. Default False.",
        ),
    ] = False,
):
    """Run a benchmark as specified in the yaml."""

    if not local:
        # TODO Check how snakemake storage decorators work, do we have caching locally or just remote?
        # TODO Implement remote execution using remote url from benchmark definition
        typer.echo(
            f"Error: Remote execution is not supported yet. Workflows can only be run in local mode.",
            err=True,
        )
        raise typer.Exit(code=1)

    else:
        benchmark = validate_benchmark(benchmark)

        if update and not dry:
            update_prompt = typer.confirm(
                "Are you sure you want to re-run the entire workflow?", abort=True
            )
            if not update_prompt:
                raise typer.Abort()

        # Controlling resource allocation with Snakemake is tricky
        # -c only controls the number of parallelism for the Snakemake scheduler
        # bioinfo methods are not designed with limited resources in mind (most)
        # Future: Create yaml for communicating resources for individual methods
        typer.echo("Running benchmark...")
        success = workflow.run_workflow(
            benchmark, cores=cores, update=update, dryrun=dry
        )

        if success:
            typer.echo(
                "Benchmark run has finished successfully.", color=typer.colors.GREEN
            )
        else:
            typer.echo("Benchmark run has failed.", err=True, color=typer.colors.RED)

        raise typer.Exit(code=0 if success else 1)


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
    input: Annotated[
        Optional[str],
        typer.Option(
            "--input",
            "-i",
            help="Run on inputs from input directory.",
        ),
    ] = None,
    example: Annotated[
        Optional[bool],
        typer.Option(
            "--example",
            "-e",
            help="Run on remote example inputs only.",
        ),
    ] = None,
    all: Annotated[
        Optional[bool],
        typer.Option(
            "--all",
            "-a",
            help="Run on all valid benchmark inputs.",
        ),
    ] = None,
    update: Annotated[
        bool,
        typer.Option(
            "--update",
            "-u",
            help="Force re-run execution for all modules and stages.",
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
):
    """
    Run a specific module on various datasets or custom inputs. This command does not have a default behavior. You must explicitly choose one of the following options:

    1. `--input`: Provide a path to a custom directory to use as the input dataset.
    2. `--example`: Set this flag to execute the module on a remote example dataset.
    3. `--all`: Set this flag to run the module on all available remote datasets.

    Note: You must select one of these options for the command to run.
    """

    behaviours = {"input": input, "example": example, "all": all}

    non_none_behaviours = {
        key: value for key, value in behaviours.items() if value is not None
    }
    if len(non_none_behaviours) == 0:
        typer.echo(
            "Error: At least one option must be specified. Use '--input', '--example', or '--all'.",
            err=True,
        )
        raise typer.Exit(code=1)
    elif len(non_none_behaviours) >= 2:
        typer.echo(
            "Error: Only one of '--input', '--example', or '--all' should be set. Please choose only one option.",
            err=True,
        )
        raise typer.Exit(code=1)
    else:
        # Construct a message specifying which option is set
        behaviour = list(non_none_behaviours)[0]

        if behaviour == "example" or behaviour == "all":
            if behaviour == "example":
                typer.echo("Running module on a predefined remote example dataset.")
            if behaviour == "all":
                typer.echo("Running module on all available remote datasets.")

            # TODO Check how snakemake storage decorators work, do we have caching locally or just remote?
            # TODO Implement remote execution using remote url from benchmark definition
            typer.echo(
                "Error: Remote execution is not supported yet. Workflows can only be run in local mode.",
                err=True,
            )
            raise typer.Exit(code=1)
        else:
            typer.echo(f"Running module on a dataset provided in a custom directory.")
            benchmark = validate_benchmark(benchmark)

            if update and not dry:
                update_prompt = typer.confirm(
                    "Are you sure you want to re-run the entire workflow?", abort=True
                )
                if not update_prompt:
                    raise typer.Abort()

            benchmark_nodes = benchmark.get_nodes_by_module_id(module_id=module)
            if len(benchmark_nodes) > 0:
                typer.echo(
                    f"Found {len(benchmark_nodes)} workflow nodes for module {module}."
                )
                typer.echo("Running module benchmark...")

                # Check if input path exists and is a directory
                if os.path.exists(input) and os.path.isdir(input):
                    benchmark_datasets = benchmark.get_benchmark_datasets()

                    # Check available files in input to figure out what dataset are we processing
                    # if we're given the initial dataset module to process, then we know
                    if module in benchmark_datasets:
                        dataset = module

                    # else we try to figure the dataset based on the files present in the input directory
                    else:
                        files = os.listdir(input)
                        base_names = [file.split(".")[0] for file in files]
                        dataset = next(
                            (d for d in benchmark_datasets if d in base_names), None
                        )

                    if dataset is not None:
                        # Check if input directory contains all necessary input files
                        required_inputs = list(
                            map(lambda node: node.get_inputs(), benchmark_nodes)
                        )
                        required_inputs = list(chain.from_iterable(required_inputs))
                        required_input_files = list(
                            set([os.path.basename(path) for path in required_inputs])
                        )
                        required_input_files = [
                            file.format(dataset=dataset)
                            for file in required_input_files
                        ]

                        input_files = os.listdir(input)
                        missing_files = [
                            file
                            for file in required_input_files
                            if file not in input_files
                        ]

                        if len(missing_files) == 0:
                            for benchmark_node in benchmark_nodes:
                                # When running a single module, it doesn't have sense to make parallelism level (cores) configurable
                                success = workflow.run_node_workflow(
                                    node=benchmark_node,
                                    input_dir=Path(input),
                                    dataset=dataset,
                                    cores=1,
                                    update=update,
                                    dryrun=dry,
                                )

                                if success:
                                    typer.echo(
                                        "Module run has finished successfully.",
                                        color=typer.colors.GREEN,
                                    )
                                else:
                                    typer.echo(
                                        "Module run has failed.",
                                        err=True,
                                        color=typer.colors.RED,
                                    )

                                raise typer.Exit(code=0 if success else 1)

                        else:
                            typer.echo(
                                f"Error: The following required input files are missing from the input directory: {missing_files}.",
                                err=True,
                                color=typer.colors.RED,
                            )

                            raise typer.Exit(code=1)

                    else:
                        typer.echo(
                            f"Error: Could not infer the appropriate dataset to run the node workflow on based on the files available in `{input}`. None of the available datasets {benchmark_datasets} match the base names of the files.",
                            err=True,
                            color=typer.colors.RED,
                        )

                        raise typer.Exit(code=1)
                else:
                    typer.echo(
                        f"Error: Input directory does not exist or is not a valid directory: `{input}`",
                        err=True,
                        color=typer.colors.RED,
                    )

                    raise typer.Exit(code=1)

            else:
                typer.echo(
                    f"Error: Could not find module with id `{module}` in benchmark definition",
                    err=True,
                    color=typer.colors.RED,
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
                f"Error: Failed to parse YAML as a valid OmniBenchmark: {str(e)}.",
                err=True,
                color=typer.colors.RED,
            )
            raise typer.Exit(code=1)

        except yaml.YAMLError as e:
            typer.echo(
                f"Error: YAML file format error: {e}.", err=True, color=typer.colors.RED
            )
            raise typer.Exit(code=1)

        except FileNotFoundError:
            typer.echo(
                "Error: Benchmark YAML file not found.",
                err=True,
                color=typer.colors.RED,
            )
            raise typer.Exit(code=1)

        except Exception as e:
            typer.echo(
                f"Error: An unexpected error occurred: {e}",
                err=True,
                color=typer.colors.RED,
            )
            raise typer.Exit(code=1)

    else:
        typer.echo(
            "Error: Invalid benchmark input. Please provide a valid YAML file path.",
            err=True,
            color=typer.colors.RED,
        )
        raise typer.Exit(code=1)
