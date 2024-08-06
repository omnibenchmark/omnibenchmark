""" General utils functions"""

from linkml_runtime.loaders import yaml_loader
from pathlib import Path
from typing import List, Union, Any


def as_list(input: Union[List, Any]):
    return input if isinstance(input, List) else [input]


def parse_instance(path: Path, target_class):
    """Load a model of target_class from a file."""
    return yaml_loader.load(str(path), target_class)
