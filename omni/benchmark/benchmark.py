"""Code to extract benchmark infos and configuration from benchmark yaml file"""

from pathlib import Path
from typing import Mapping, List

import omni_schema.datamodel.omni_schema as model
import yaml

from ..utils import parse_instance


def _dict_from_yaml(bench_yaml: Path) -> Mapping:
    with open(bench_yaml, "r") as file:
        return yaml.safe_load(file)


def load_benchmark_from_yaml(bench_yaml: Path) -> model.Benchmark:
    return parse_instance(bench_yaml, model.Benchmark)


def get_steps(benchmark: model.Benchmark) -> List[str]:
    return [step.id for step in benchmark.steps]


def get_step_by_id(benchmark: model.Benchmark, step_id: str) -> model.Step:
    step = list(filter(lambda step: step.id == step_id, benchmark.steps))
    if len(step) < 1:
        raise ValueError(
            f"No step with id {step_id} found. Avaliable ids are: {get_steps(benchmark)}"
        )
    return step[0]


def get_input_collection(benchmark: model.Benchmark, step_id: str) -> List:
    step = get_step_by_id(benchmark, step_id)
    input_ids = step.inputs.entries


def get_explicit_inputs(
    benchmark: model.Benchmark, step_id: str, test: bool = True
) -> Mapping:
    pass


def get_explicit_outputs(benchmark: model.Benchmark, step_id: str) -> Mapping:
    pass


def get_available_parameter(benchmark: model.Benchmark, step_id: str) -> Mapping:
    pass


def _get_benchmark_outputs(benchmark: model.Benchmark) -> List[model.IOFile]:
    return {out_id: out_file for out_id, out_file in benchmark.steps}
