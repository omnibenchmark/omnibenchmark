"""cli commands related to benchmark/module execution and start"""

import os
import sys
from itertools import chain
from pathlib import Path

import click


@click.group(name="run")
@click.pass_context
def run(ctx):
    """Run benchmarks or benchmark modules."""
    ctx.ensure_object(dict)


@run.command(name="benchmark")
@click.pass_context
@click.option(
    "-b",
    "--benchmark",
    help="Path to benchmark yaml file or benchmark id.",
    required=True,
    envvar="OB_BENCHMARK",
    type=click.Path(exists=True),
)
@click.option(
    "-p",
    "--threads",
    help="The parallelism level for the workflow scheduler.",
    type=int,
    default=1,
)
@click.option(
    "-u",
    "--update",
    help="Force re-run execution for all modules and stages.",
    is_flag=True,
    default=False,
)
@click.option("-d", "--dry", help="Dry run.", is_flag=True, default=False)
@click.option(
    "-l",
    "--local",
    help="Execute and store results locally. Default False.",
    is_flag=True,
    default=False,
)
def run_benchmark(ctx, benchmark, threads, update, dry, local):
    """Run a benchmark as specified in the yaml."""
    ctx.ensure_object(dict)
    from omni.benchmark import Benchmark
    from omni.io.utils import remote_storage_snakemake_args
    from omni.workflow.snakemake import SnakemakeEngine
    from omni.workflow.workflow import WorkflowEngine

    benchmark = validate_benchmark(benchmark)
    workflow: WorkflowEngine = SnakemakeEngine()

    if update and not dry:
        update_prompt = click.confirm(
            "Are you sure you want to re-run the entire workflow?", abort=True
        )
        if not update_prompt:
            raise click.Abort()

    if not local:
        storage_options = remote_storage_snakemake_args(benchmark)
    else:
        storage_options = {}

    # Controlling resource allocation with Snakemake is tricky
    # -c only controls the number of parallelism for the Snakemake scheduler
    # bioinfo methods are not designed with limited resources in mind (most)
    # Future: Create yaml for communicating resources for individual methods
    click.echo("Running benchmark...")
    success = workflow.run_workflow(
        benchmark, cores=threads, update=update, dryrun=dry, **storage_options
    )

    if success:
        click.echo("Benchmark run has finished successfully.")
    else:
        click.echo("Benchmark run has failed.", err=True)

    # raise click.Exit(code=0 if success else 1)


@run.command(name="module")
@click.pass_context
@click.option(
    "-b",
    "--benchmark",
    help="Path to benchmark yaml file or benchmark id.",
    required=True,
    envvar="OB_BENCHMARK",
    type=click.Path(exists=True),
)
@click.option("-m", "--module", help="Module id to execute", type=str, required=True)
@click.option(
    "-i",
    "--input_dir",
    help="Path to the folder with the appropriate input files.",
    type=click.Path(exists=True, writable=True),
    default=None,
)
@click.option("-d", "--dry", help="Dry run.", is_flag=True, default=False)
@click.option(
    "-u",
    "--update",
    help="Force re-run execution for all modules and stages.",
    is_flag=True,
    default=False,
)
def run_module(ctx, benchmark, module, input_dir, dry, update):
    """
    Run a specific module that is part of the benchmark.
    """
    behaviours = {"input": input_dir, "example": None, "all": None}

    non_none_behaviours = {
        key: value for key, value in behaviours.items() if value is not None
    }
    if len(non_none_behaviours) >= 2:
        click.echo(
            "Error: Only one of '--input_dir', '--example', or '--all' should be set. Please choose only one option.",
            err=True,
        )
        sys.exit(1)  # raise click.Exit(code=1)
    else:
        # Construct a message specifying which option is set
        behaviour = list(non_none_behaviours)[0] if non_none_behaviours else None

        if behaviour == "example" or behaviour == "all":
            if behaviour == "example":
                click.echo("Running module on a predefined remote example dataset.")
            if behaviour == "all":
                click.echo("Running module on all available remote datasets.")

            # TODO Check how snakemake storage decorators work, do we have caching locally or just remote?
            # TODO Implement remote execution using remote url from benchmark definition
            click.echo(
                "Error: Remote execution is not supported yet. Workflows can only be run in local mode.",
                err=True,
            )
            sys.exit(1)  # raise click.Exit(code=1)
        else:
            click.echo(f"Running module on a local dataset.")
            benchmark = validate_benchmark(benchmark)

            benchmark_nodes = benchmark.get_nodes_by_module_id(module_id=module)
            is_entrypoint_module = all(
                [node.is_entrypoint() for node in benchmark_nodes]
            )
            if len(benchmark_nodes) > 0:
                click.echo(
                    f"Found {len(benchmark_nodes)} workflow nodes for module {module}."
                )

                if not is_entrypoint_module and len(non_none_behaviours) == 0:
                    click.echo(
                        "Error: At least one option must be specified. Use '--input_dir', '--example', or '--all'.",
                        err=True,
                    )
                    sys.exit(1)  # raise click.Exit(code=1)
                elif input_dir is None:
                    input_dir = os.getcwd()

                click.echo("Running module benchmark...")

                # Check if input path exists and is a directory
                if os.path.exists(input_dir) and os.path.isdir(input_dir):
                    benchmark_datasets = benchmark.get_benchmark_datasets()

                    # Check available files in input to figure out what dataset are we processing
                    # if we're given the initial dataset module to process, then we know
                    if module in benchmark_datasets:
                        dataset = module

                    # else we try to figure the dataset based on the files present in the input directory
                    else:
                        files = os.listdir(input_dir)
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

                        input_files = os.listdir(input_dir)
                        missing_files = [
                            file
                            for file in required_input_files
                            if file not in input_files
                        ]

                        if len(missing_files) == 0:
                            from omni.workflow.snakemake import SnakemakeEngine
                            from omni.workflow.workflow import WorkflowEngine

                            workflow: WorkflowEngine = SnakemakeEngine()
                            for benchmark_node in benchmark_nodes:
                                # When running a single module, it doesn't have sense to make parallelism level (cores) configurable
                                success = workflow.run_node_workflow(
                                    node=benchmark_node,
                                    input_dir=Path(input_dir),
                                    dataset=dataset,
                                    cores=1,
                                    update=update,
                                    dryrun=dry,
                                )

                                if success:
                                    click.echo(
                                        "Module run has finished successfully.",
                                    )
                                else:
                                    click.echo(
                                        "Module run has failed.",
                                        err=True,
                                    )

                                sys.exit(
                                    0 if success else 1
                                )  # raise click.Exit(code=0 if success else 1)

                        else:
                            click.echo(
                                f"Error: The following required input files are missing from the input directory: {missing_files}.",
                                err=True,
                            )

                            sys.exit(1)  # raise click.Exit(code=1)

                    else:
                        click.echo(
                            f"Error: Could not infer the appropriate dataset to run the node workflow on based on the files available in `{input_dir}`. None of the available datasets {benchmark_datasets} match the base names of the files.",
                            err=True,
                        )

                        sys.exit(1)  # raise click.Exit(code=1)
                else:
                    click.echo(
                        f"Error: Input directory does not exist or is not a valid directory: `{input_dir}`",
                        err=True,
                    )

                    sys.exit(1)  # raise click.Exit(code=1)

            else:
                click.echo(
                    f"Error: Could not find module with id `{module}` in benchmark definition",
                    err=True,
                )
                sys.exit(1)  # raise click.Exit(code=1)


@run.command(no_args_is_help=True, name="validate")
@click.pass_context
@click.option(
    "-b",
    "--benchmark",
    help="Path to benchmark yaml file or benchmark id.",
    envvar="OB_BENCHMARK",
    type=click.Path(exists=True),
)
def validate_yaml(ctx, benchmark):
    """Validate a benchmark yaml."""
    click.echo("Validating a benchmark yaml.")
    benchmark = validate_benchmark(benchmark)


## to validate the YAML
def validate_benchmark(benchmark_file: str):
    from pathlib import Path

    import yaml

    from omni.benchmark import Benchmark

    if benchmark_file.endswith(".yaml") or benchmark_file.endswith(".yml"):
        try:
            with open(benchmark_file, "r") as file:
                yaml.safe_load(file)
                benchmark = Benchmark(Path(benchmark_file))
                click.echo("Benchmark YAML file integrity check passed.")

                return benchmark

        except ValueError as e:
            click.echo(
                f"Error: Failed to parse YAML as a valid OmniBenchmark: {str(e)}.",
                err=True,
            )
            sys.exit(1)  # raise click.Exit(code=1)

        except yaml.YAMLError as e:
            click.echo(f"Error: YAML file format error: {e}.", err=True)
            sys.exit(1)  # raise click.Exit(code=1)

        except FileNotFoundError:
            click.echo(
                "Error: Benchmark YAML file not found.",
                err=True,
            )
            sys.exit(1)  # raise click.Exit(code=1)

        except Exception as e:
            click.echo(
                f"Error: An unexpected error occurred: {e}",
                err=True,
            )
            sys.exit(1)  # raise click.Exit(code=1)

    else:
        click.echo(
            "Error: Invalid benchmark input. Please provide a valid YAML file path.",
            err=True,
        )
        sys.exit(1)  # raise click.Exit(code=1)
