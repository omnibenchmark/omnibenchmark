import re
from itertools import takewhile
from pathlib import Path
from typing import List, Set, Tuple, Union, NamedTuple, Dict

from omnibenchmark.model import MetricCollector

from omnibenchmark.benchmark import BenchmarkNode, Benchmark
from omnibenchmark.utils import format_mc_output


class Wildcards(NamedTuple):
    pre: str
    post: str
    dataset: str


def format_performance_file(node: BenchmarkNode, out_dir: str) -> str:
    """Provides a benchmark performance path for a node"""

    is_initial = node.is_initial()
    if is_initial:
        return out_dir + "/{stage}/{module}/{params}/{dataset}_performance.txt"
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
    benchmark: Benchmark, collector: MetricCollector, return_as_dict: bool = False
) -> Union[Dict[str, str], List[str]]:
    """Formats collector inputs that will be expanded according to Snakemake's engine"""
    # Handle both string and IOFile inputs
    # We should probably decide on one or the other, keeping both is messy.
    implicit_inputs = []
    for i in collector.inputs:
        if isinstance(i, str):
            implicit_inputs.append(i)
        else:
            implicit_inputs.append(i.id)
    explicit_inputs = benchmark.get_explicit_input(implicit_inputs)

    if not return_as_dict:
        inputs = list(explicit_inputs.values())
        return inputs
    else:
        return explicit_inputs

    # print(f'Collector Input: {collector.name}: {inputs}')


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
    benchmark: Benchmark, wildcards: Wildcards, return_as_dict: bool = False
) -> Union[Dict[str, str], List[str]]:
    """Formats benchmark inputs that will be expanded according to Snakemake's engine"""

    pre = wildcards.pre
    post = wildcards.post
    dataset = wildcards.dataset

    assert not dataset.__contains__(
        "/"
    ), "Dataset wildcard should not contain a path prefix"

    nodes = benchmark.get_nodes()
    stage_ids = set(benchmark.get_stage_ids())

    pre_stages = _extract_stages_from_path(pre, stage_ids)
    matched = _match_node_format(pre_stages[-1]) if len(pre_stages) > 0 else []

    if len(matched) != 3:
        return {} if return_as_dict else []

    after_stage_id, _, _ = matched
    stage_id, module_id, param_id = _match_node_format(post)

    node_hash = hash(BenchmarkNode.to_id(stage_id, module_id, param_id, after_stage_id))
    matching_node = next((node for node in nodes if hash(node) == node_hash), None)
    if not matching_node:
        return {} if return_as_dict else []

    node_inputs = matching_node.get_inputs_dict()

    inputs = _match_inputs(node_inputs, pre_stages, pre, dataset)

    if return_as_dict:
        return inputs
    else:
        return list(inputs.values())


def _extract_stages_from_path(
    path: str, known_stage_ids: Set[str]
) -> List[Tuple[str, ...]]:
    def is_known_stage(part: str) -> bool:
        return part in known_stage_ids

    def find_sub_parts(start_index: int) -> Tuple[str, ...]:
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


def _match_node_format(to_match: Union[str, Tuple[str, ...]]) -> Tuple[str, str, str]:
    if isinstance(to_match, str):
        # TODO(daninci): Fix type checking
        parts = [part for part in Path(to_match).parts]
        to_match = tuple(parts)

    stage_id = to_match[0]
    module_id = to_match[1]
    param_id = to_match[2]

    return stage_id, module_id, param_id


def _match_input_module(input: str, stages: List[Tuple[str, ...]], dataset: str) -> str:
    expected_input_module = input.split("{pre}/")[1].split("/{module}")[0]
    matching_stage = next(
        (tup for tup in stages if tup[0] == expected_input_module), None
    )
    if not matching_stage:
        raise RuntimeError(f"Could not find matching stage for {input} in {stages}")

    if len(matching_stage) < 3:
        raise RuntimeError(f"{stages} has wrong length")

    assert len(matching_stage) >= 3
    matched_module = matching_stage[1]

    input = input.replace("{module}", matched_module)
    input = input.replace("{dataset}", dataset)
    if "{params}" in input:
        # Use first parameter value if available, otherwise default to "default"
        matched_params = matching_stage[2] if len(matching_stage) > 2 else "default"
        input = input.replace("{params}", matched_params)

    return input


def _match_input_prefix(input: str, pre: str) -> str:
    stage = f"/{input.split('/')[1]}"
    matched_prefix = pre.split(stage)[0]
    formatted_input = input.format(pre=matched_prefix)

    # print(f'Input: {input} Prefix: {pre} Formatted: {formatted_input}')
    return formatted_input


def _match_inputs(
    inputs: Dict[str, str], stages: List[Tuple[str, ...]], pre: str, dataset: str
) -> Dict[str, str]:
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
