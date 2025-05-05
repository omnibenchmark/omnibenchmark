import re
from itertools import takewhile
from pathlib import Path
from typing import List, Set, Tuple, Union, NamedTuple, Dict

from omni_schema.datamodel.omni_schema import MetricCollector

from omnibenchmark.benchmark import BenchmarkNode, Benchmark
from omnibenchmark.utils import format_mc_output


class Wildcards(NamedTuple):
    pre: str
    post: str
    dataset: str


def format_performance_file(node: BenchmarkNode) -> str:
    """Provides a benchmark performance path for a node"""

    is_initial = node.is_initial()
    if is_initial:
        return "out/{stage}/{module}/{params}/{dataset}_performance.txt"
    else:
        return "{pre}/{post}/{dataset}_performance.txt"


def format_output_templates_to_be_expanded(node: BenchmarkNode) -> List[str]:
    """Formats node outputs that will be expanded according to Snakemake's engine"""
    outputs = node.get_outputs()
    is_initial = node.is_initial()

    if not is_initial:
        outputs = [re.sub(r"\/.*\/", "/{post}/", o, count=1) for o in outputs]

    # print(f'Output: {node.stage_id}: {outputs}')
    return outputs


def format_metric_collector_input(
    benchmark: Benchmark, collector: MetricCollector, return_as_dict=False
) -> Dict[str, str] | List[str]:
    """Formats collector inputs that will be expanded according to Snakemake's engine"""
    implicit_inputs = [i.id for i in collector.inputs]
    explicit_inputs = benchmark.get_explicit_input(implicit_inputs)

    if not return_as_dict:
        inputs = list(explicit_inputs.values())
    else:
        inputs = explicit_inputs

    # print(f'Collector Input: {collector.name}: {inputs}')
    return inputs


def format_metric_collector_output(
    out_dir: Path, collector: MetricCollector
) -> List[str]:
    """Formats collector outputs that will be expanded according to Snakemake's engine"""
    outputs = []
    for o in collector.outputs:
        outputs.append(format_mc_output(o, out_dir, collector.id))

    # print(f'Collector Output: {collector.name}: {outputs}')
    return outputs


def format_input_templates_to_be_expanded(
    benchmark: Benchmark, wildcards: Wildcards, return_as_dict=False
) -> Dict[str, str] | List[str]:
    """Formats benchmark inputs that will be expanded according to Snakemake's engine"""

    pre = wildcards.pre
    post = wildcards.post
    dataset = wildcards.dataset

    nodes = benchmark.get_nodes()
    stage_ids = benchmark.get_stage_ids()

    pre_stages = _extract_stages_from_path(pre, stage_ids)
    after_stage_id, _, _ = (
        _match_node_format(pre_stages[-1]) if len(pre_stages) > 0 else None
    )

    stage_id, module_id, param_id = _match_node_format(post)

    node_hash = hash(BenchmarkNode.to_id(stage_id, module_id, param_id, after_stage_id))
    matching_node = next((node for node in nodes if hash(node) == node_hash), None)
    if matching_node:
        node_inputs = matching_node.get_inputs_dict()

        inputs = _match_inputs(node_inputs, pre_stages, pre, dataset)

        # print(f'Inputs: {stage_id} {module_id} {param_id}: {inputs}')
        if return_as_dict:
            return inputs
        else:
            return inputs.values()
    else:
        return {} if return_as_dict else []


def _extract_stages_from_path(
    path: str, known_stage_ids: Set[str]
) -> List[Union[str, tuple]]:
    def is_known_stage(part: str) -> bool:
        return part in known_stage_ids

    def find_sub_parts(start_index: int) -> Union[str, tuple]:
        sub_parts = list(
            takewhile(lambda x: x not in known_stage_ids, parts[start_index + 1 :])
        )
        return tuple(parts[start_index : start_index + 1 + len(sub_parts)])

    if not path:
        return []

    parts = [part for part in Path(path).parts]
    stages = [
        find_sub_parts(part_idx)
        for part_idx, part in enumerate(parts)
        if is_known_stage(part)
    ]

    return stages


def _match_node_format(to_match: Union[str, tuple]) -> Tuple[str, str, str]:
    if isinstance(to_match, str):
        # TODO(daninci): Fix type checking
        to_match = [part for part in Path(to_match).parts]

    stage_id = to_match[0]
    module_id = to_match[1]
    param_id = to_match[2]

    return stage_id, module_id, param_id


def _match_input_module(input: str, stages: List[Tuple[str]], dataset: str) -> str:
    expected_input_module = input.split("{pre}/")[1].split("/{module}")[0]
    matching_stage = next(
        (tup for tup in stages if tup[0] == expected_input_module), None
    )

    if matching_stage:
        matched_module = matching_stage[1]

        input = input.replace("{module}", matched_module)
        input = input.replace("{dataset}", dataset)
        if "{params}" in input:
            matched_params = next(
                (x for x in matching_stage[2:] if "param" or "default" in x), None
            )
            input = input.replace("{params}", matched_params)

        return input
    else:
        raise RuntimeError(f"Could not find matching stage for {input} in {stages}")


def _match_input_prefix(input: str, pre: str) -> str:
    stage = f"/{input.split('/')[1]}"
    matched_prefix = pre.split(stage)[0]
    formatted_input = input.format(pre=matched_prefix)

    # print(f'Input: {input} Prefix: {pre} Formatted: {formatted_input}')
    return formatted_input


def _match_inputs(
    inputs: dict[str, str], stages: List[Tuple[str]], pre: str, dataset: str
) -> dict[str, str]:
    all_matched = True

    formatted_inputs = {}
    for key, input in inputs.items():
        formatted_input = _match_input_module(input, stages, dataset)
        if not formatted_input:
            all_matched = False
            break
        else:
            formatted_input = _match_input_prefix(formatted_input, pre)
            formatted_inputs[key] = formatted_input

    return formatted_inputs if all_matched else {}
