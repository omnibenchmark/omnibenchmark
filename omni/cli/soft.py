"""cli commands related to software management"""

from typing_extensions import Annotated

from omni.software import easybuild_backend
from omni.software import conda_backend
from omni.software import common

import typer

import os

# import logging

cli = typer.Typer(add_completion = True,  no_args_is_help = True, pretty_exceptions_short = True)

# @cli.callback()
# def main(verbose: bool = typer.Option(False, "--verbose", "-v")):
#     app = typer.Typer(pretty_exceptions_short = True)
#     if verbose:
#         app = typer.Typer(pretty_exceptions_short = False)

# @cli.callback()
# def main(pretty: bool = typer.Option(False, "--pretty", "-p")):
#     app = typer.Typer(pretty_exceptions_short = False)
#     if pretty:
#         app = typer.Typer(pretty_exceptions_enable = True)

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
    typer.echo(f"Install software for {easyconfig} within a Singularity container.")
    
    ## check lmod
    easybuild_backend.export_lmod_env_vars()
    # easybuild_backend.check_module_tool()

    ## check easybuild
    easybuild_backend.check_easybuild_status()
    
    ## do
    envmodule_name = easybuild_backend.get_envmodule_name_from_easyconfig(easyconfig, workdir = os.getcwd())    
    easybuild_backend.create_definition_file(
        easyconfig = easyconfig,
        singularity_recipe = 'Singularity_' + easyconfig + '.txt',
        envmodule = envmodule_name, nthreads = str(len(os.sched_getaffinity(0))),
        templates_path = os.path.dirname(easybuild_backend.__file__))
    print('skipping the shell part now')
    

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


@cli.command("check")
def check_singularity(
    singularity: Annotated[
        str,
        typer.Option(
            "--sif",
            "-s",
            help="Path to the singularity SIF.",
        ),
    ],
    envmodule_name: Annotated[
        str,
        typer.Option(
            "--envmodulename",
            "-e",
            help="Envmodule name.",
        ),
    ],
):
    """Check whether the envmodule {envmodule} is loadable (installed) in a given singularity sif {sif}."""
    typer.echo(f"Checks whether the envmodule {envmodule} is loadable and installed in a given singularity sif {sif}.")
