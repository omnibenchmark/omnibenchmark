"""cli commands related to software management"""

from pathlib import Path
from typing_extensions import Annotated

from omni.software import easybuild_backend as eb
from omni.software import conda_backend
from omni.software import common

import typer

import os, sys

# import logging

cli = typer.Typer(
    add_completion=True,
    no_args_is_help=True,
    pretty_exceptions_short=False,
    rich_markup_mode=None,
    pretty_exceptions_enable=False,
    help="Manage benchmark-specific software installations.",
)

sing_cli = typer.Typer(
    add_completion=True, no_args_is_help=True, pretty_exceptions_short=True
)
cli.add_typer(
    sing_cli,
    name="singularity",
    help="Manage singularity- (apptainer-) based software installations. Uses easybuild to build software.",
)

conda_cli = typer.Typer(
    add_completion=True, no_args_is_help=True, pretty_exceptions_short=True
)
cli.add_typer(
    conda_cli,
    name="conda",
    help="Manage conda-based software installations. Does not use easybuild.",
)

module_cli = typer.Typer(
    add_completion=True, no_args_is_help=True, pretty_exceptions_short=True
)
cli.add_typer(
    module_cli,
    name="envmodules",
    help="Manage envmodules-based software installations. Uses easybuild to build software.",
)

docker_cli = typer.Typer(
    add_completion=True, no_args_is_help=True, pretty_exceptions_short=True
)
cli.add_typer(
    docker_cli,
    name="docker",
    help="Manage docker based software installations. Uses easybuild to build software.",
)

## singularity #####################################################################################


@sing_cli.command("build")
def singularity_build(
    easyconfig: Annotated[
        str,
        typer.Option(
            "--easyconfig",
            "-e",
            help="Easyconfig",
        ),
    ],
):
    """Build a singularity (fakeroot) image for a given easyconfig."""
    typer.echo(
        f"Installing software for {easyconfig} within a Singularity container. It will take some time."
    )

    if common.check_easybuild_status().returncode != 0:
        raise ("ERROR: Easybuild not installed")
    if common.check_singularity_status().returncode != 0:
        raise ("ERROR: Singularity not installed")

    ## check the easyconfig exists
    try:
        fp = eb.get_easyconfig_full_path(easyconfig=easyconfig)
    except:
        typer.echo("ERROR: easyconfig not found.\n", err=True, color=typer.colors.RED)
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

    eb.singularity_build(singularity_recipe=singularity_recipe, easyconfig=easyconfig)
    typer.echo("DONE: recipe and image built for " + singularity_recipe)


@sing_cli.command("push")
def singularity_push(
    docker_username: Annotated[
        str,
        typer.Option(
            "--docker_username",
            "-u",
            help="Docker username",
        ),
    ],
    docker_password: Annotated[
        str,
        typer.Option(
            "--docker_password",
            "-p",
            help="Docker password (token)",
        ),
    ],
    sif: Annotated[
        Path,
        typer.Option(
            "--sif",
            "-s",
            help="Path to the Singularity SIF file",
        ),
    ],
    oras: Annotated[
        str,
        typer.Option(
            "--oras",
            "-o",
            help="Registry's ORAS static URL, for instance oras://registry.mygitlab.ch/myuser/myproject:mytag",
        ),
    ],
):
    """Pushes a singularity SIF file to an ORAS-compatible registry."""
    typer.echo(f"Pushing {sif} to the registry {oras}.")

    eb.push_to_registry(
        sif=sif,
        docker_username=docker_username,
        docker_password=docker_password,
        oras=oras,
    )
    typer.echo("DONE\n.")


## envmodules ########################################################################################


@module_cli.command("build")
def envmodules_build(
    easyconfig: Annotated[
        str,
        typer.Option(
            "--easyconfig",
            "-e",
            help="Easyconfig",
        ),
    ],
    threads: Annotated[
        int,
        typer.Option("--threads", "-p", help="Number of threads"),
    ] = 2,
):
    """Build a given easyconfig (and generates the relevant envmodules)."""
    typer.echo(
        f"Installing software for {easyconfig} using easybuild. It will take some time."
    )

    if common.check_easybuild_status().returncode != 0:
        raise ("ERROR: Easybuild not installed")

    ## check the easyconfig exists
    try:
        fp = eb.get_easyconfig_full_path(easyconfig=easyconfig)
    except:
        print("ERROR: easyconfig not found.\n")
        sys.exit()

    eb.easybuild_easyconfig(easyconfig=easyconfig, threads=threads)
    typer.echo("DONE")


# @module_cli.command('list')
# def envmodules_list():
#     """Lists available modules"""
#     typer.echo(f"Listing modules")

#     eb.list_available_modules("*")

#     print('DONE')

# @module_cli.command('load')
# def envmodule_load(
#     module: Annotated[
#         str,
#         typer.Option(
#             "--module",
#             "-m",
#             help="Module name",
#         ),
#     ]
# ):
#     """Loads an envmodule"""
#     typer.echo(f"Loading module {module}")

#     if common.check_easybuild_status().returncode != 0:
#         raise('ERROR: Easybuild not installed')

#     eb.load_module_api(mod_name = module)
#     # print(module)

#     print('DONE')


## conda #############################################################################################


@conda_cli.command("pin")
def pin_conda_env(
    conda_env: Annotated[
        str,
        typer.Option(
            "--env",
            "-e",
            help="Path to the conda env file.",
        ),
    ],
):
    """Pin all conda env-related dependencies versions using snakedeploy."""
    typer.echo(f"Pinning {conda_env} via snakedeploy. It will take some time.")
    conda_backend.pin_conda_envs(conda_env)
    typer.echo(f"\nDONE: Pinned {conda_env}\n")


## general stuff ######################################################################################


@cli.command("check")
def check(
    what: Annotated[
        str,
        typer.Option(
            "--what",
            "-w",
            help="""Binary/functionality to check: \n
               --what singularity : singularity \n
               --what module      : module tool, typically lmod \n 
               --what easybuild   : easybuild \n
               --what conda       : conda \n
               --what docker      : docker \n""",
        ),
    ],
):
    """Check whether the component {what} is available."""
    typer.echo(
        f"Checking software stack handlers / backends (singularity, easybuild, etc)."
    )

    if what == "easybuild":
        ret = common.check_easybuild_status()
    elif what == "module":
        # eb.export_lmod_env_vars()
        ret = common.check_lmod_status()
    elif what == "singularity":
        ret = common.check_singularity_status()
    elif what == "conda":
        ret = common.check_conda_status()
    elif what == "docker":
        ret = common.check_docker_status()
    else:
        raise typer.BadParameter(
            "Bad `--what` value. Please check help (`ob software check --help`)."
        )
    if ret.returncode == 0:
        typer.echo("OK: " + ret.stdout)
    else:
        typer.echo("Failed: " + ret.returncode)


## docker


@docker_cli.command("build")
def docker_build(
    easyconfig: Annotated[
        str,
        typer.Option(
            "--easyconfig",
            "-e",
            help="Easyconfig",
        ),
    ],
):
    """Build a docker image for a given easyconfig."""
    typer.echo(
        f"Installing software for {easyconfig} within a docker container. It will take some time."
    )

    if common.check_easybuild_status().returncode != 0:
        raise ("ERROR: Easybuild not installed")
    if common.check_docker_status().returncode != 0:
        raise ("ERROR: Docker not installed")

    ## check the easyconfig exists
    try:
        fp = eb.get_easyconfig_full_path(easyconfig=easyconfig)
    except:
        typer.echo("ERROR: easyconfig not found.\n", err=True, color=typer.colors.RED)
        sys.exit()

    docker_recipe = "Dockerfile_" + easyconfig + ".txt"
    envmodule_name = eb.get_envmodule_name_from_easyconfig(easyconfig)
    eb.create_dockerfile(
        easyconfig=easyconfig,
        dockerfile=docker_recipe,
        envmodule=envmodule_name,
        nthreads=str(len(os.sched_getaffinity(0))),
    )
    typer.echo("DONE: dockerfile built for " + docker_recipe)
    eb.docker_build(
        dockerfile=docker_recipe, easyconfig=easyconfig, name=envmodule_name
    )
    typer.echo("DONE: docker image built for " + docker_recipe)
