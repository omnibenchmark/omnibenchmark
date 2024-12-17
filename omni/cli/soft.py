"""cli commands related to software management"""

import os
import sys
from pathlib import Path

import click
import yaml

from omni.cli.utils.logging import debug_option, logger


# import logging
@click.group()
@click.pass_context
@debug_option
def software(ctx):
    """Manage and install benchmark-specific software."""
    ctx.ensure_object(dict)


## singularity #####################################################################################


@click.group
@click.pass_context
def singularity(ctx):
    """Manage and install software using Singularity."""
    ctx.ensure_object(dict)


software.add_command(singularity)


@singularity.command(name="build", no_args_is_help=True)
@click.option(
    "-e",
    "--easyconfig",
    help="Easyconfig.",
    required=True,
)
def singularity_build(easyconfig):
    """Build a singularity (fakeroot) image for a given easyconfig."""
    logger.info(
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
        logger.error("ERROR: easyconfig not found.")
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

    logger.info(
        "DONE: singularity recipe written for "
        + singularity_recipe
        + "\nDOING: building the image"
    )

    sb = eb.singularity_build(
        singularity_recipe=singularity_recipe, easyconfig=easyconfig
    )
    if sb.returncode != 0:
        logger.error("ERROR: " + sb.stderr)
    if sb.returncode == 0:
        logger.info(sb.stdout)
        logger.info("DONE: singularity image built for " + singularity_recipe)


@singularity.command(name="prepare")
@click.option(
    "-b",
    "--benchmark",
    help="Benchmark YAML.",
    type=click.Path(exists=True),
    required=True,
    envvar="OB_BENCHMARK",
)
def singularity_prepare(benchmark):
    """Build all singularity (fakeroot) images needed for a benchmark."""
    logger.info(
        f"Installing software for {benchmark} using Singularity containers. It will take some time."
    )
    from omni.benchmark import Benchmark
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
            logger.error("ERROR: easyconfig not found.")
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

        logger.info(
            "DONE: singularity recipe written for "
            + singularity_recipe
            + "\nDOING: building the image"
        )

        sb = eb.singularity_build(
            singularity_recipe=singularity_recipe, easyconfig=easyconfig
        )
        if sb.returncode != 0:
            logger.error("ERROR: " + sb.stderr)
        if sb.returncode == 0:
            logger.info(sb.stdout)
            logger.info("DONE: singularity image built for " + singularity_recipe)
        logger.info("DONE: singularity images built.")


@singularity.command(name="push", no_args_is_help=True)
@click.option(
    "-u", "--docker_username", help="Docker username.", type=str, required=True
)
@click.option(
    "-p", "--docker_password", help="Docker password.", type=str, required=True
)
@click.option(
    "-s",
    "--sif",
    help="Path to the Singularity SIF file.",
    type=click.Path(writable=True),
    required=True,
)
@click.option(
    "-o",
    "--oras",
    help="Registry's ORAS static URL, for instance oras://registry.mygitlab.ch/myuser/myproject:mytag.",
    type=str,
    required=True,
)
def singularity_push(docker_username, docker_password, sif, oras):
    """Pushes a singularity SIF file to an ORAS-compatible registry."""
    logger.info(f"Pushing {sif} to the registry {oras}.")

    from omni.software import easybuild_backend as eb

    eb.push_to_registry(
        sif=sif,
        docker_username=docker_username,
        docker_password=docker_password,
        oras=oras,
    )
    logger.info("DONE\n.")


## envmodules ########################################################################################


@click.group
@click.pass_context
def module(ctx):
    """Manage and install software using Easybuild"""
    ctx.ensure_object(dict)


software.add_command(module)


@module.command(name="build", no_args_is_help=True)
@click.option(
    "-e",
    "--easyconfig",
    help="Easyconfig.",
    required=True,
)
@click.option("-p", "--threads", type=int, default=2, help="Number of threads.")
def envmodules_build(easyconfig, threads):
    """Build a given easyconfig (and generates the relevant envmodules)."""
    logger.info(
        f"\nInstalling software for {easyconfig} using easybuild. It will take some time.\n"
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
        logger.error("ERROR: " + p.stderr)
    if p.returncode == 0:
        logger.info(p.stdout)
    logger.info("DONE: built " + easyconfig)


@module.command(name="prepare", no_args_is_help=True)
@click.option(
    "-b",
    "--benchmark",
    help="Benchmark YAML.",
    type=click.Path(exists=True),
    required=True,
    envvar="OB_BENCHMARK",
)
@click.option("-p", "--threads", default=2, help="Number of threads.", type=int)
def envmodules_prepare(benchmark, threads):
    """Build all envmodules needed for a given benchmark YAML."""
    logger.info(
        f"Installing software for {benchmark} using envmodules. It will take some time."
    )
    from omni.benchmark import Benchmark
    from omni.software import common
    from omni.software import easybuild_backend as eb

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
            logger.error(
                "ERROR: easyconfig not found.\n",
            )
            sys.exit()

        p = eb.easybuild_easyconfig(easyconfig=easyconfig, threads=threads)
        if p.returncode != 0:
            logger.error("ERROR: " + p.stderr)
        if p.returncode == 0:
            logger.info(p.stdout)
            logger.info("DONE: built " + easyconfig)
        logger.info("DONE: built all easyconfigs")


## conda #############################################################################################


@click.group
@click.pass_context
def conda(ctx):
    """Manage and install software using conda."""
    ctx.ensure_object(dict)


software.add_command(conda)


@conda.command(name="pin")
@click.option(
    "-e",
    "--env",
    "conda_env",
    help="Conda env YAML.",
    type=click.Path(writable=True),
    required=True,
)
def pin_conda_env(conda_env):
    """Pin all conda env-related dependencies versions using snakedeploy."""
    logger.info(f"Pinning {conda_env} via snakedeploy. It will take some time.")
    from omni.software import common, conda_backend

    conda_backend.pin_conda_envs(conda_env)
    logger.info(f"\nDONE: Pinned {conda_env}\n")


@conda.command(name="prepare")
@click.option(
    "-b",
    "--benchmark",
    help="Benchmark YAML.",
    type=click.Path(exists=True),
    required=True,
    envvar="OB_BENCHMARK",
)
def conda_prepare(benchmark):
    """Pin all conda envs needed for a given benchmark YAML."""
    logger.info(f"Pinning conda envs for {benchmark}. It will take some time.")
    from omni.benchmark import Benchmark
    from omni.software import common, conda_backend

    if common.check_conda_status().returncode != 0:
        raise ("ERROR: conda not installed")

    with open(benchmark, "r") as fh:
        yaml.safe_load(fh)
        benchmark = Benchmark(Path(benchmark))

    for conda in benchmark.get_conda_envs():
        if not os.path.isfile(os.path.join("envs", conda)):
            logger.error(
                "ERROR: theconda env file at "
                + os.path.join("envs", conda)
                + "does not exist."
            )
            sys.exit()
        conda_backend.pin_conda_envs(os.path.join("envs", conda))
    logger.info("DONE: pinned all conda envs.")


## general stuff ######################################################################################


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
    logger.info(
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
        logger.info("OK: " + ret.stdout)
    else:
        logger.error("Failed: " + ret.stdout + ret.stderr)


software.add_command(check)
