""" click-based CLI """

import click

# from pathlib import Path

# from omni.benchmark import Benchmark
# from omni.software import easybuild_backend as eb
# from omni.software import conda_backend
# from omni.software import common
# import yaml

import os, sys


@click.group()
@click.option("--debug/--no-debug", default=False)
@click.pass_context
def cli(ctx, debug):
    ctx.ensure_object(dict)

    ctx.obj["DEBUG"] = debug


@click.group()
@click.pass_context
def info(ctx):
    ctx.ensure_object(dict)


@click.group()
@click.pass_context
def run(ctx):
    ctx.ensure_object(dict)


@click.group()
@click.pass_context
def software(ctx):
    ctx.ensure_object(dict)


@click.group()
@click.pass_context
def storage(ctx):
    ctx.ensure_object(dict)


cli.add_command(info)
cli.add_command(run)
cli.add_command(software)
cli.add_command(storage)

## info start        ###########################################################################

## info end          ###########################################################################


## run start         ###########################################################################


@run.command(no_args_is_help=True, name="benchmark")
@click.pass_context
@click.option("-b", "--benchmark", help="Path to benchmark yaml file or benchmark id.")
@click.option(
    "-p", "--threads", help="The parallelism level for the workflow scheduler."
)
@click.option(
    "-u", "--update", help="Force re-run execution for all modules and stages."
)
@click.option("-d", "--dry", help="Dry run.")
def run_benchmark(ctx, benchmark, threads, update, dry):
    """Run a benchmark as specified in the yaml."""
    ctx.ensure_object(dict)
    from omni.benchmark import Benchmark
    from omni.io.utils import remote_storage_snakemake_args
    from omni.workflow.snakemake import SnakemakeEngine
    from omni.workflow.workflow import WorkflowEngine

    benchmark = validate_benchmark(benchmark)
    workflow: WorkflowEngine = SnakemakeEngine()

    # if update and not dry:
    #     update_prompt = click.confirm(
    #         "Are you sure you want to re-run the entire workflow?", abort=True
    #     )
    #     if not update_prompt:
    #         raise click.Abort()

    # if not local:
    #     storage_options = remote_storage_snakemake_args(benchmark)
    # else:
    #     storage_options = {}

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


@run.command(no_args_is_help=True, name="module")
@click.pass_context
@click.option("-b", "--benchmark", help="Path to benchmark yaml file or benchmark id.")
@click.option("-m", "--module", help="Module id to execute")
@click.option(
    "-i", "--input_dir", help="Path to the folder with the appropriate input files."
)
@click.option("-d", "--dry", help="Dry run.")
def run_module(ctx, benchmark, module, input, dry):
    """
    Run a specific module on inputs present at a given folder.
    """
    example = False
    all = False
    behaviours = {"input": input, "example": example, "all": all}

    non_none_behaviours = {
        key: value for key, value in behaviours.items() if value is not None
    }
    if len(non_none_behaviours) == 0:
        click.echo(
            "Error: At least one option must be specified. Use '--input', '--example', or '--all'.",
            err=True,
        )
        sys.exit(1)  # raise click.Exit(code=1)

    elif len(non_none_behaviours) >= 2:
        click.echo(
            "Error: Only one of '--input', '--example', or '--all' should be set. Please choose only one option.",
            err=True,
        )
        sys.exit(1)  # raise click.Exit(code=1)
    else:
        # Construct a message specifying which option is set
        behaviour = list(non_none_behaviours)[0]

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
            click.echo(f"Running module on a dataset provided in a custom directory.")
            benchmark = validate_benchmark(benchmark)

            # if update and not dry:
            #     update_prompt = click.confirm(
            #         "Are you sure you want to re-run the entire workflow?", abort=True
            #     )
            #     if not update_prompt:
            #         raise click.Abort()

            benchmark_nodes = benchmark.get_nodes_by_module_id(module_id=module)
            if len(benchmark_nodes) > 0:
                click.echo(
                    f"Found {len(benchmark_nodes)} workflow nodes for module {module}."
                )
                click.echo("Running module benchmark...")

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
                            f"Error: Could not infer the appropriate dataset to run the node workflow on based on the files available in `{input}`. None of the available datasets {benchmark_datasets} match the base names of the files.",
                            err=True,
                        )

                        sys.exit(1)  # raise click.Exit(code=1)
                else:
                    click.echo(
                        f"Error: Input directory does not exist or is not a valid directory: `{input}`",
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
@click.option("-b", "--benchmark", help="Path to benchmark yaml file or benchmark id.")
def validate_yaml(ctx, benchmark):
    """Validate a benchmark yaml."""
    click.echo("Validating a benchmark yaml.")
    benchmark = validate_benchmark(benchmark)


## to validate the YAML
def validate_benchmark(benchmark_file: str):
    from omni.benchmark import Benchmark
    from pathlib import Path
    import yaml

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


## run end           ###########################################################################


## software start    ###########################################################################


@click.group
@click.pass_context
def singularity(ctx):
    ctx.ensure_object(dict)


software.add_command(singularity)


@singularity.command(name="build", no_args_is_help=True)
@click.option("-e", "--easyconfig", help="Easyconfig.")
def singularity_build(easyconfig):
    """Build a singularity (fakeroot) image for a given easyconfig."""
    click.echo(
        f"Installing software for {easyconfig} within a Singularity container. It will take some time."
    )
    from omni.software import common
    from omni.software import easybuild_backend as eb

    if common.check_easybuild_status().returncode != 0:
        raise ("ERROR: Easybuild not installed")
    if common.check_singularity_status().returncode != 0:
        raise ("ERROR: Singularity not installed")

    ## check the easyconfig exists
    try:
        fp = eb.get_easyconfig_full_path(easyconfig=easyconfig)
    except:
        click.echo("ERROR: easyconfig not found.\n", err=True)
        sys.exit()

    ## do
    singularity_recipe = "Singularity_" + easyconfig + ".txt"
    envmodule_name = eb.get_envmodule_name_from_easyconfig(easyconfig)
    eb.create_definition_file(
        easyconfig=easyconfig,
        singularity_recipe=singularity_recipe,
        envmodule=envmodule_name,
        nthreads=str(len(os.sched_getaffinity(0))),
    )

    click.echo(
        "DONE: singularity recipe written for "
        + singularity_recipe
        + "\nDOING: building the image"
    )

    sb = eb.singularity_build(
        singularity_recipe=singularity_recipe, easyconfig=easyconfig
    )
    if sb.returncode != 0:
        click.echo("ERROR: " + sb.stderr)
    if sb.returncode == 0:
        click.echo(sb.stdout)
        click.echo("DONE: singularity image built for " + singularity_recipe)


@singularity.command(name="prepare", no_args_is_help=True)
@click.option("-b", "--benchmark", help="Benchmark YAML.")
def singularity_prepare(benchmark):
    """Build all singularity (fakeroot) images needed for a benchmark."""
    click.echo(
        f"Installing software for {benchmark} using Singularity containers. It will take some time."
    )
    from omni.software import common
    from omni.software import easybuild_backend as eb

    if common.check_easybuild_status().returncode != 0:
        raise ("ERROR: Easybuild not installed")
    if common.check_singularity_status().returncode != 0:
        raise ("ERROR: Singularity not installed")

    with open(benchmark, "r") as fh:
        yaml.safe_load(fh)
        benchmark = Benchmark(Path(benchmark))

    for easyconfig in benchmark.get_easyconfigs():
        print(easyconfig)
        ## check the easyconfig exists
        try:
            fp = eb.get_easyconfig_full_path(easyconfig=easyconfig)
        except:
            click.echo("ERROR: easyconfig not found.\n", err=True)
            sys.exit()

        ## do
        singularity_recipe = "Singularity_" + easyconfig + ".txt"
        envmodule_name = eb.get_envmodule_name_from_easyconfig(easyconfig)
        eb.create_definition_file(
            easyconfig=easyconfig,
            singularity_recipe=singularity_recipe,
            envmodule=envmodule_name,
            nthreads=str(len(os.sched_getaffinity(0))),
        )

        click.echo(
            "DONE: singularity recipe written for "
            + singularity_recipe
            + "\nDOING: building the image"
        )

        sb = eb.singularity_build(
            singularity_recipe=singularity_recipe, easyconfig=easyconfig
        )
        if sb.returncode != 0:
            click.echo("ERROR: " + sb.stderr)
        if sb.returncode == 0:
            click.echo(sb.stdout)
            click.echo("DONE: singularity image built for " + singularity_recipe)
        click.echo("DONE: singularity images built.")


@singularity.command(name="push", no_args_is_help=True)
@click.option("-u", "--docker_username", help="Docker username.")
@click.option("-p", "--docker_password", help="Docker password.")
@click.option("-s", "--sif", help="Path to the Singularity SIF file.")
@click.option(
    "-o",
    "--oras",
    help="Registry's ORAS static URL, for instance oras://registry.mygitlab.ch/myuser/myproject:mytag.",
)
def singularity_push(docker_username, docker_password, sif, oras):
    """Pushes a singularity SIF file to an ORAS-compatible registry."""
    click.echo(f"Pushing {sif} to the registry {oras}.")

    from omni.software import easybuild_backend as eb

    eb.push_to_registry(
        sif=sif,
        docker_username=docker_username,
        docker_password=docker_password,
        oras=oras,
    )
    click.echo("DONE\n.")


@click.group
@click.pass_context
def module(ctx):
    ctx.ensure_object(dict)


software.add_command(module)


@module.command(name="build", no_args_is_help=True)
@click.option("-e", "--easyconfig", help="Easyconfig.")
@click.option("-p", "--threads", default=2, help="Number of threads.")
def envmodules_build(easyconfig):
    """Build a given easyconfig (and generates the relevant envmodules)."""
    click.echo(
        f"Installing software for {easyconfig} using easybuild. It will take some time."
    )

    from omni.software import common
    from omni.software import easybuild_backend as eb

    if common.check_easybuild_status().returncode != 0:
        raise ("ERROR: Easybuild not installed")

    ## check the easyconfig exists
    try:
        fp = eb.get_easyconfig_full_path(easyconfig=easyconfig)
    except:
        print("ERROR: easyconfig not found.\n")
        sys.exit()

    p = eb.easybuild_easyconfig(easyconfig=easyconfig, threads=threads)
    if p.returncode != 0:
        click.echo("ERROR: " + p.stderr)
    if p.returncode == 0:
        click.echo(p.stdout)
    click.echo("DONE: built " + easyconfig)


@module.command(name="prepare", no_args_is_help=True)
@click.option("-b", "--benchmark", help="Benchmark YAML.")
@click.option("-p", "--threads", default=2, help="Number of threads.")
def envmodules_prepare(benchmark, threads):
    """Build all envmodules needed for a given benchmark YAML."""
    click.echo(
        f"Installing software for {benchmark} using envmodules. It will take some time."
    )

    if common.check_easybuild_status().returncode != 0:
        raise ("ERROR: Easybuild not installed")
    if common.check_lmod_status().returncode != 0:
        raise ("ERROR: lmod not installed")

    with open(benchmark, "r") as fh:
        yaml.safe_load(fh)
        benchmark = Benchmark(Path(benchmark))

    for easyconfig in benchmark.get_easyconfigs():
        print(easyconfig)
        ## check the easyconfig exists
        try:
            fp = eb.get_easyconfig_full_path(easyconfig=easyconfig)
        except:
            click.echo(
                "ERROR: easyconfig not found.\n",
                err=True,
            )
            sys.exit()

        p = eb.easybuild_easyconfig(easyconfig=easyconfig, threads=threads)
        if p.returncode != 0:
            click.echo("ERROR: " + p.stderr)
            if p.returncode == 0:
                click.echo(p.stdout)
            click.echo("DONE: built " + easyconfig)
        click.echo("DONE: built all easyconfigs")


@click.group
@click.pass_context
def conda(ctx):
    ctx.ensure_object(dict)


software.add_command(conda)


@module.command(name="pin", no_args_is_help=True)
@click.option("-e", "--env", help="Conda env YAML.")
def pin_conda_env(conda_env):
    """Pin all conda env-related dependencies versions using snakedeploy."""
    click.echo(f"Pinning {conda_env} via snakedeploy. It will take some time.")
    from omni.software import common
    from omni.software import conda_backend

    conda_backend.pin_conda_envs(conda_env)
    click.echo(f"\nDONE: Pinned {conda_env}\n")


@module.command(name="prepare", no_args_is_help=True)
@click.option("-b", "--benchmark", help="Benchmark YAML.")
def conda_prepare(benchmark):
    """Pin all conda envs needed for a given benchmark YAML."""
    click.echo(f"Pinning conda envs for {benchmark}. It will take some time.")
    from omni.software import common
    from omni.software import conda_backend

    if common.check_conda_status().returncode != 0:
        raise ("ERROR: conda not installed")

    with open(benchmark, "r") as fh:
        yaml.safe_load(fh)
        benchmark = Benchmark(Path(benchmark))

    for conda in benchmark.get_conda_envs():
        if not os.path.isfile(os.path.join("envs", conda)):
            click.echo(
                "ERROR: theconda env file at "
                + os.path.join("envs", conda)
                + "does not exist."
            )
            sys.exit()
        conda_backend.pin_conda_envs(os.path.join("envs", conda))
    click.echo("DONE: pinned all conda envs.")


@click.command(no_args_is_help=True)
@click.pass_context
@click.option(
    "-w",
    "--what",
    help="""Binary/functionality to check: \n
               --what singularity : singularity \n
               --what module      : module tool, typically lmod \n 
               --what easybuild   : easybuild \n
               --what conda       : conda \n""",
)
def check(ctx, what):
    """Check whether the component {what} is available."""
    click.echo(
        f"Checking software stack handlers / backends (singularity, easybuild, etc)."
    )
    from omni.software import common

    if what == "easybuild":
        ret = common.check_easybuild_status()
    elif what == "module":
        # eb.export_lmod_env_vars()
        ret = common.check_lmod_status()
    elif what == "singularity":
        ret = common.check_singularity_status()
    elif what == "conda":
        ret = common.check_conda_status()
    # elif what == "docker":
    #     ret = common.check_docker_status()
    else:
        raise click.BadParameter(
            "Bad `--what` value. Please check help (`ob software check --help`)."
        )
    if ret.returncode == 0:
        click.echo("OK: " + ret.stdout)
    else:
        click.echo("Failed: " + ret.stdout + ret.stderr)


software.add_command(check)

## software end      ###########################################################################

## storage start     ###########################################################################

## storage end       ###########################################################################

if __name__ == "__main__":
    cli(obj={})
