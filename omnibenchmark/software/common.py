"""
Common omniblock software management logics, regardless of conda/easybuild/system.
Calls workflow capabilities to reuse schema / yaml parsing logics
"""

import subprocess


def check_call(command: str, use_shell: bool = True) -> subprocess.CompletedProcess:
    return subprocess.run(
        command,
        text=True,
        capture_output=True,
        check=False,
        shell=use_shell,
    )


def check_lmod_status() -> subprocess.CompletedProcess:
    ret = subprocess.run(
        ["type", "module"], text=True, capture_output=True, check=False, shell=True
    )
    return ret


def check_singularity_status() -> subprocess.CompletedProcess:
    return check_call("singularity --version")


def check_easybuild_status() -> subprocess.CompletedProcess:
    return check_call("eb --version")


def check_conda_status() -> subprocess.CompletedProcess:
    return check_call("conda --version")
