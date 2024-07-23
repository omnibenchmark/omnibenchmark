""" General utils functions"""
from linkml_runtime.loaders import yaml_loader
from pathlib import Path
from typing import List, Union, Any
import re


def as_list(input: Union[List, Any]):
    return input if isinstance(input, List) else [input]


def parse_instance(path: Path, target_class):
    """Load a model of target_class from a file."""
    return yaml_loader.load(str(path), target_class)


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
