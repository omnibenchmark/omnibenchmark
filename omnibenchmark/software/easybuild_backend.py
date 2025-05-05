#!/usr/bin/env python

"""
Easybuild-powered software management, mostly by omniblock. Also includes singularity image generation.

Izaskun Mallona
Started 5th June 2024
"""

import subprocess
import os
import os.path as op

from easybuild.tools.module_naming_scheme.utilities import det_full_ec_version
from easybuild.framework.easyconfig.tools import det_easyconfig_paths, parse_easyconfigs
from easybuild.tools.options import set_up_configuration


from importlib import resources as impresources
from . import templates

HOME = op.expanduser("~")
## careful, 'all' will be prepended unless exported, shell-wise, the MODULEPATH/all path as MODULEPATH
##   when running omb
# MODULEPATH=op.join(HOME, 'micromamba', 'envs', 'omnibenchmark', 'modulefiles', 'Core')
MODULEPATH = op.join(HOME, ".local", "easybuild", "modules")
ROBOTPATH = op.join(
    HOME, "micromamba", "envs", "omnibenchmark", "easybuild", "easyconfigs"
)

os.makedirs(MODULEPATH, exist_ok=True)

## shell-based stuff, partly to be replaced by direct eb API calls -------------------------------------


def generate_default_easybuild_config_arguments(
    modulepath: str = MODULEPATH,
    sourcepath: str = ROBOTPATH + ":" + os.getcwd(),
) -> str:
    args = "--installpath-modules=" + modulepath + " --sourcepath=" + sourcepath

    return args


set_up_configuration(
    args=["--info"] + generate_default_easybuild_config_arguments().split(" "),
    silent=True,
)


def construct_easybuild_easyconfig_command(easyconfig: str, threads: int = 2) -> str:
    """
    Constructs an easybuild command to build an easyconfig

     Args:
     - easyconfig (str): the easyconfig name. Doesn't have to be a full path. But readable from the robots path.
     - threads (int): number of threads to build the software
    """

    args = generate_default_easybuild_config_arguments()
    cmd = (
        """eb %(easyconfig)s --robot --parallel=%(threads)s %(args)s --detect-loaded-modules=unload --check-ebroot-env-vars=unset"""
        % {"easyconfig": op.basename(easyconfig), "threads": threads, "args": args}
    )

    return cmd


def easybuild_easyconfig(
    easyconfig: str,
    threads: int = 2,
) -> subprocess.CompletedProcess:
    """
    Easybuilds an easyconfig

     Args:
     - easyconfig (str): the easyconfig name. Doesn't have to be a full path. But readable from the robots path.
     - threads (int): number of threads to build the software
    """

    cmd = construct_easybuild_easyconfig_command(
        easyconfig=op.basename(easyconfig), threads=threads
    )

    # try:
    ret = subprocess.run(
        cmd.split(" "), shell=False, capture_output=True, check=False, text=True
    )
    # print(ret.stdout)
    # print(ret.stderr)
    return ret
    # except subprocess.CalledProcessError as exc:
    #     return ("ERROR easybuild failed:", exc.returncode, exc.output)


def parse_easyconfig(ec_fn: str) -> list:
    """
    Find and parse an easyconfig with specified filename,
    and return parsed easyconfig file (an EasyConfig instance).

    Args:
    - ec_fn (str): easyconfig filename. Doesn't have to be a full path. But readable from the robots path

    """
    # opts, _ = set_up_configuration(args = generate_default_easybuild_config_arguments(workdir).split(),
    #                                silent = True)

    ec_path = det_easyconfig_paths([ec_fn])[0]

    # parse easyconfig file;
    # the 'parse_easyconfigs' function expects a list of tuples,
    # where the second item indicates whether or not the easyconfig file was automatically generated or not
    ec_dicts, _ = parse_easyconfigs([(ec_path, False)])

    # only retain first parsed easyconfig
    return ec_path, ec_dicts[0]["ec"]


def get_easyconfig_full_path(easyconfig: str) -> str:
    """
    Finds the easyconfig full path
    Args:
    - easyconfig (str): easyconfig filename. Doesn't have to be a full path. But readable from the robots path
    """

    try:
        ec_path, ec = parse_easyconfig(op.basename(easyconfig))
        return ec_path
    except Exception:
        raise FileNotFoundError("\nERROR: easyconfig not found.\n")


def get_envmodule_name_from_easyconfig(easyconfig: str) -> str:
    """
    Returns the (standard) envmodulename from an easyconfig file.
    Args:
    - easyconfig (str): easyconfig filename. Doesn't have to be a full path. But readable from the robots path
    """
    ec_path, ec = parse_easyconfig(op.basename(easyconfig))
    return os.path.join(ec["name"], det_full_ec_version(ec))


# def construct_easybuild_easyconfig_command(
#     easyconfig : str, workdir :str =  os.getcwd(), threads : int = 2, containerize : bool =False, container_build_image : bool =False
# ) -> str:


#     # args = generate_default_easybuild_config_arguments(workdir = workdir)
#     cmd = """eb %(easyconfig)s --robot --parallel=%(threads)s \
#               --detect-loaded-modules=unload --check-ebroot-env-vars=unset""" % {
#         "easyconfig": easyconfig,
#         # 'args' : args,
#         "threads": threads,
#     }
#     if containerize:
#         cmd += (
#             " --container-config bootstrap=localimage,from=example.sif --experimental"
#         )
#         if container_build_image:
#             cmd += " --container-build-image"
#     return cmd


def create_definition_file(
    easyconfig: str, singularity_recipe: str, envmodule: str, nthreads: int = 2
) -> None:
    """
    Materializes a Singularity recipe to build a given easyconfig, using a Singularity recipe template.
    Args:
    - easyconfig (str): easyconfig filename. Doesn't have to be a full path. But readable from the robots path
    - singularity_recipe (str): path to a Singularity template
    - envmodule (str): the envmodule name
    - nthreads (int): number of threads for easybuild
    """
    template = impresources.files(templates) / "ubuntu_jammy.txt"
    with open(template, "rt") as ubuntu, open(singularity_recipe, "w") as sing:
        for line in ubuntu.read().split("\n"):
            if "EASYCONFIG" in line:
                line = line.replace("EASYCONFIG", op.basename(easyconfig))
            if "ENVMODULENAME" in line:
                line = line.replace("ENVMODULENAME", envmodule)
            if "EASYBUILDNTHREADSINT" in line:
                line = line.replace("EASYBUILDNTHREADSINT", nthreads)
            sing.write(line + "\n")


def singularity_build(
    easyconfig: str, singularity_recipe: str
) -> subprocess.CompletedProcess:
    """
    Builds a Singularity recipe
    Args:
    - easyconfig (str): easyconfig filename. Doesn't have to be a full path. But readable from the robots path
    - singularity_recipe (str): path to a Singularity template
    """
    image_name = op.basename(easyconfig) + ".sif"

    cmd = (
        "singularity build --force --fakeroot " + image_name + " " + singularity_recipe
    )
    # output = subprocess.check_output(
    #     cmd, stderr=subprocess.STDOUT, shell=True, universal_newlines=True
    # )
    output = subprocess.run(
        cmd, shell=True, text=True, capture_output=True, check=False
    )
    return output


## untested,drafted 06 Aug 2024
def singularity_push(
    sif: str, docker_username: str, docker_password: str, oras_url: str
) -> subprocess.CompletedProcess | None:
    cmd = f"""singularity push --docker-username {docker_username} \
                 --docker-password {docker_password} \
                 {sif} \
                 {oras_url}"""
    try:
        subprocess.run(
            cmd.split(" "), shell=False, text=True, capture_output=True, check=True
        )
    except Exception as exc:
        raise RuntimeError("ERROR singularity build failed:", exc)
