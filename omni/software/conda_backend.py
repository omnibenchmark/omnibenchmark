#!/usr/bin/env python

from snakedeploy import conda as snakeconda

"""
Conda-based software management, per benchmark yaml (and not per omniblock)
"""

def pin_conda_envs(benchmark_yaml):
    """
    Pins a conda env file using snakedeploy
    """
    snakeconda.pin_conda_envs(conda_env_paths = [benchmark_yaml],
                              conda_frontend = "mamba", warn_on_error= False)

def check_node_conda(benchmark_yaml, stage_id, module_id):
    """
    Check if the env.yaml and the pin.txt for that conda env file exist locally.
    If so: valid, continue (Snakemake will install, if applicable).
    If not:
       - env.yaml not present: error
       - env.yaml present but pin.txt not present: error
    """
    return('not implemented')

# forget about those, we assume conda exists
# def install_miniconda():

# def install_mamba():
#     try:
#         subprocess.call(["mamba", "--help"])
#     except FileNotFoundError:
#         print("ERROR mamba not found; please install it system-wise")
#     try:
#         subprocess.call(["conda", "--help"])
#     except FileNotFoundError:
#         print("ERROR conda not found, please install it")
#     try:
#         cmd = 'conda install conda-forge::mamba'
#         output = subprocess.check_output(
#             cmd, stderr = subprocess.STDOUT, shell = True,
#             universal_newlines = True)
#     except subprocess.CalledProcessError as exc:
#         return("ERROR mamba install failed:", exc.returncode, exc.output)
#     else:
#         return("LOG mamba install output: \n{}\n".format(output))
