"""General utils functions"""

import os

from linkml_runtime.loaders import yaml_loader
import subprocess
from pathlib import Path
from typing import List, Union, Any
import re
import yaml
import platform

from omni_schema.datamodel.omni_schema import IOFile


def try_avail_envmodule(module_name: str) -> bool:
    env = {}
    env.update(os.environ)

    command = f"""
    . "$LMOD_PKG"/init/profile ;
    module purge ;
    module avail {module_name}"""
    result = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=True,
        text=True,
        env=env,
    )

    return "No module(s) or extension(s) found!" not in result.stderr


def as_list(input: Union[List, Any]):
    return input if isinstance(input, List) else [input]


def parse_instance(path: Path, target_class) -> Any:
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


def format_mc_output(output: IOFile, out_dir: Path, collector_id: str):
    if output.path:
        o = output.path.replace("{input}", str(out_dir))
        o = o.replace("{name}", collector_id)
        return o
    else:
        return str(out_dir / output.id)
