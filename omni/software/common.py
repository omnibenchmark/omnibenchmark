#!/usr/bin/env python

"""
Common omniblock software management logics, regardless of conda/easybuild/system. 
Calls workflow capabilities to reuse schema / yaml parsing logics
"""

# import omni_workflow import main

import subprocess


def check_call(command: str, use_shell: bool = True) -> subprocess.CompletedProcess:
    ret = subprocess.run(
        command.split(" "),
        text=True,
        capture_output=True,
        check=True,
        shell=use_shell,
    )
    return ret


def check_lmod_status() -> subprocess.CompletedProcess:
    return check_call("type module", use_shell=True)


def check_singularity_status() -> subprocess.CompletedProcess:
    return check_call("singularity --version")


def check_easybuild_status() -> subprocess.CompletedProcess:
    return check_call("eb --version")


def check_conda_status() -> subprocess.CompletedProcess:
    return check_call("conda --version")


def check_docker_status() -> subprocess.CompletedProcess:
    return check_call("docker --version")
