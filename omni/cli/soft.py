"""cli commands related to software management"""

from typing_extensions import Annotated

from omni.software import easybuild_backend
from omni.software import conda_backend
from omni.software import common

import typer

import os

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
    typer.echo(f"Installing software for {easyconfig} within a Singularity container.")
    
    ## check lmod
    easybuild_backend.export_lmod_env_vars()
    # easybuild_backend.check_module_tool()

    ## check easybuild
    easybuild_backend.check_easybuild_status()
    
    ## do
    singularity_recipe = 'Singularity_' + easyconfig + '.txt'
    envmodule_name = easybuild_backend.get_envmodule_name_from_easyconfig(easyconfig, workdir = os.getcwd())    
    easybuild_backend.create_definition_file(
        easyconfig = easyconfig,
        singularity_recipe = singularity_recipe,
        envmodule = envmodule_name, nthreads = str(len(os.sched_getaffinity(0))),
        templates_path = os.path.dirname(easybuild_backend.__file__))

    easybuild_backend.singularity_build(easyconfig = easyconfig, singularity_recipe = singularity_recipe)
    

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
    typer.echo(f"Pinning {conda_env} via snakedeploy.")
    conda_backend.pin_conda_envs(conda_env)


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
        easybuild_backend.check_easybuild_status()
    elif what == 'lmod':
        easybuild_backend.export_lmod_env_vars()
        easybuild_backend.check_lmod_status()
    elif what == 'singularity':
        easybuild_backend.check_singularity_status()
    else:
        return('not implemented')
