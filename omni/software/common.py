#!/usr/bin/env python

"""
Common omniblock software management logics, regardless of conda/easybuild/system. 
Calls workflow capabilities to reuse schema / yaml parsing logics
"""

# import omni_workflow import main

import subprocess


def check_call(command: str, use_shell: bool = False) -> subprocess.CompletedProcess:
    try:
        ret = subprocess.run(
            command.split(" "),
            text=True,
            capture_output=True,
            check=True,
            shell=use_shell,
        )
        return ret
    except Exception as e:
        return print("ERROR:", e)


def check_lmod_status() -> subprocess.CompletedProcess:
    return check_call("type module", use_shell=True)


def check_singularity_status() -> subprocess.CompletedProcess:
    return check_call("singularity --version")


def check_easybuild_status() -> subprocess.CompletedProcess:
    return check_call("eb --version")


def check_conda_status() -> subprocess.CompletedProcess:
    return check_call("conda --version")
