"""Code to extract benchmark infos and configuration from benchmark yaml file"""
from pathlib import Path
import yaml
from typing import Mapping, List

def _dict_from_yaml(bench_yaml: Path) -> Mapping:
    with open(bench_yaml, 'r') as file:
        return yaml.safe_load(file)

def get_stages(bench_dict: Mapping) -> List :
    return [stage for stage in bench_dict.get("steps").keys()]

def get_modules(bench_dict: Mapping) -> List:
    stages = get_stages(bench_dict)
    