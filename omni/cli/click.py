""" click-based CLI """

import click

# from pathlib import Path

# from omni.benchmark import Benchmark
# from omni.software import easybuild_backend as eb
# from omni.software import conda_backend
# from omni.software import common
# import yaml

# import os, sys


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
        click.echo("ERROR: easyconfig not found.\n", err=True, color=click.colors.RED)
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
            click.echo(
                "ERROR: easyconfig not found.\n", err=True, color=click.colors.RED
            )
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
                "ERROR: easyconfig not found.\n", err=True, color=click.colors.RED
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
