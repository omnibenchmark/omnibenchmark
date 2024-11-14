""" General utils functions"""
import logging
import os

from linkml_runtime.loaders import yaml_loader
import subprocess
from pathlib import Path
from typing import List, Union, Any, Optional
import re
import yaml
import platform
# from env_modules_python import module
from subprocess import PIPE, Popen
import sys 


## stolen from https://github.com/TACC/Lmod/blob/main/init/env_modules_python.py.in
def module(command, *arguments, **kwargs):
    """
    Execute a regular Lmod command and apply environment changes to
    the current Python environment (i.e. os.environ).
    
    In case len(arguments) == 1 the string will be split on whitespace into 
    separate arguments. Pass a list of strings to avoid this.
    
    Raises an exception in case Lmod execution returned a non-zero
    exit code.
    
    Use with keyword argument show_environ_updates=True to show the actual
    changes made to os.environ (mostly for debugging).
    
    Examples:
    module('list')
    module('load', 'gcc')
    module('load', 'gcc cmake')
    module('load', 'gcc cmake', show_environ_updates=True)
    """
    numArgs = len(arguments)
    A = ['@PKG@/libexec/lmod', 'python', command]
    if (numArgs == 1):
        A += arguments[0].split()
    else:
        A += list(arguments)

    proc           = Popen(A, stdout=PIPE, stderr=PIPE)
    status         = proc.returncode 
    stdout, stderr = proc.communicate()
    err_out        = sys.stderr
    if (os.environ.get('LMOD_REDIRECT','@redirect@') != 'no'):
        err_out=sys.stdout

    print(stderr.decode(),file=err_out)
    
    if ('show_environ_updates' in kwargs):
        print(stdout.decode())
    exec(stdout.decode())
    return status, stderr.decode()


def is_module_available() -> bool:
    try:
        # # Run `module --version` to check if `module` command is available
        # command = """
        #     source "$LMOD_PKG"/init/profile
        #     export PYTHONPATH=${PYTHONPATH}:$LMOD_DIR/../init
        #     module --version
        # """
        # result = subprocess.run(
        #     command,
        #     stdout=subprocess.PIPE,
        #     stderr=subprocess.PIPE,
        #     shell=True,
        #     text=True,
        # )

        # # If the command executes without an error, `module` is available
        # return result.returncode == 0
        
        status = module('list')
        return status[1] == 'No modules loaded\n'
    
    except Exception as e:
        logging.error(f"ERROR: Module command does not exist: {e}")
        return False


# def get_available_modules(modulepath: str = os.environ['MODULEPATH'], module_name: Optional[str] = None) -> List[str]:
#     if not is_module_available():
#         return []

#     try:
#         # # Run `module spider` to get the list of all available modules
#         # # Combine all commands into a single shell execution
#         # command = f"""
#         #     source "$LMOD_PKG"/init/profile
#         #     #export PYTHONPATH=${{PYTHONPATH}}:$LMOD_DIR/../init
#         #     export MODULEPATH={modulepath}
#         #     module spider {module_name}"""

#         # result = subprocess.run(
#         #     command,
#         #     stdout=subprocess.PIPE,
#         #     stderr=subprocess.PIPE,
#         #     shell=True,
#         #     text=True,
#         # )

#         # print(f"\n==============================")
#         # print(f"Command: {command}")
#         # print(f"Stdout: {result.stdout.strip()}")
#         # print(f"Stderr: {result.stderr.strip()}")
#         # print(f"MODULEPATH: {os.environ['MODULEPATH']}")
#         # print("==============================")

#         # Parse the output to get module names
#         # Each module is usually listed with "  <module_name>:" pattern
#         # modules = set(
#         #     re.findall(r"^\s{2}(\S+):", result.stdout + result.stderr, re.MULTILINE)
#         # )
#         # return list(sorted(modules))
#         module('spider', f"{module_name}")
#         return f"{module_name}"
    
#     except Exception as e:
#         logging.error(f"ERROR: Failed to retrieve available lmod modules: {e}")
#         return []


def try_load_envmodule(module_name: str) -> bool:
    # if not is_module_available():
    #     return False

    try:
        # command = f"""
        #     source "$LMOD_PKG"/init/profile
        #     export MODULEPATH={modulepath}
        #     # export PYTHONPATH=${{PYTHONPATH}}:$LMOD_DIR/../init
        
        #     module load {module_name}"""
        # result = subprocess.run(
        #     command,
        #     stdout=subprocess.PIPE,
        #     stderr=subprocess.PIPE,
        #     shell=True,
        #     text=True,
        # )

        # # If the command executes without an error, `module` was loaded successfully
        # return result.returncode == 0
        status = module('load', f"{module_name}")
        return True
    except:
        return not 'error' in status[1]
    


def as_list(input: Union[List, Any]):
    return input if isinstance(input, List) else [input]


def parse_instance(path: Path, target_class):
    """Load a model of target_class from a file."""

    # Due to a yaml.CLoader issue for Yaml files on Windows, https://github.com/yaml/pyyaml/issues/293
    # We will have different logic for loading the model based on the OS.
    # Basically the Windows user will have a slower loading time, because the yaml.Loader does not have the issue
    if platform.system() == "Windows":
        with path.open("r") as file:
            benchmark_yaml = yaml.load(file, yaml.SafeLoader)
            benchmark = yaml_loader.load(benchmark_yaml, target_class)

            return benchmark
    else:
        benchmark = yaml_loader.load(str(path), target_class)
        return benchmark


def merge_dict_list(list_of_dicts):
    """Merge a list of dictionaries into a single dictionary."""
    merged_dict = {
        key: value for d in list_of_dicts if d is not None for key, value in d.items()
    }

    return merged_dict


def format_name(path, prefix):
    pattern = rf"{prefix}/.+?/([^/]+)/.+?$"
    dataset = re.match(pattern, path)[1]
    new_path = path.format(dataset=dataset)

    return new_path
