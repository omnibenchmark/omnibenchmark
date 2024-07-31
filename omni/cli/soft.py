"""cli commands related to software management"""

from typing_extensions import Annotated

from omni.software import easybuild_backend as eb
from omni.software import conda_backend
from omni.software import common

import typer

import os, sys

# import logging

cli = typer.Typer(add_completion = True,  no_args_is_help = True, pretty_exceptions_short = True)

@cli.command("build")
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
    typer.echo(f"Installing software for {easyconfig} within a Singularity container. It will take some time.")
    
    ## check lmod
    eb.export_lmod_env_vars()
    # eb.check_module_tool()

    ## check easybuild
    eb.check_easybuild_status()

    ## check the easyconfig exists
    try:
        fp = eb.get_easyconfig_full_path(easyconfig = easyconfig, workdir = os.getcwd())
    except:
        print('ERROR: easyconfig not found.\n')
        sys.exit()
        
    ## do
    singularity_recipe = 'Singularity_' + easyconfig + '.txt'
    envmodule_name = eb.get_envmodule_name_from_easyconfig(easyconfig, workdir = os.getcwd())    
    eb.create_definition_file(
        easyconfig = easyconfig,
        singularity_recipe = singularity_recipe,
        envmodule = envmodule_name, nthreads = str(len(os.sched_getaffinity(0))))

    eb.singularity_build(easyconfig = easyconfig, singularity_recipe = singularity_recipe)
    print('DONE: recipe and image built for ' + singularity_recipe)

@cli.command("pin")
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
    """Pins all conda env-related dependencies versions using snakedeploy."""
    typer.echo(f"Pinning {conda_env} via snakedeploy. It will take some time.")
    conda_backend.pin_conda_envs(conda_env)
    typer.echo(f'DONE: Pinned {conda_env}\n')


## these belong to the validator, not here

@cli.command("check")
def check(
    what: Annotated[
        str,
        typer.Option(
            "--what",
            "-w",
            help="""Binary/functionality to check: 
               --what singularity : singularity binary 
               --what lmod        : module tool 
               --what easybuild   : easybuild binary""",
        ),
    ],        
):
    """Check whether the component {what} is available."""
    typer.echo(f"Checking whether the component {what} is available.")

    if what == 'easybuild':
        eb.check_easybuild_status()
    elif what == 'lmod':
        eb.export_lmod_env_vars()
        eb.check_lmod_status()
    elif what == 'singularity':
        eb.check_singularity_status()
    else:
        return('not implemented')
