from linkml_runtime.loaders import yaml_loader
from pathlib import Path


def parse_instance(path: Path, target_class):
    """Load a model of target_class from a file."""
    return yaml_loader.load(str(path), target_class)
