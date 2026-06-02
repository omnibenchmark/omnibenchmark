from pathlib import Path
from typing import List, Tuple, Optional

from omnibenchmark.dag import (
    DiGraph,
    topological_sort,
    all_simple_paths,
    NetworkXUnfeasible,
)
from omnibenchmark.model import Benchmark, Stage, ValidationError

from ._node import BenchmarkNode


def expand_stage_nodes(
    model: Benchmark,
    stage: "Stage",
    out_dir: Path,
    stage_ordering: List[str],
) -> List[BenchmarkNode]:
    nodes = []

    input_dirname = str(out_dir) if model.is_initial(stage) else "{pre}"
    stage_outputs = model.get_stage_outputs(stage).values()
    outputs = [x.replace("{input}", input_dirname) for x in stage_outputs]

    inputs_for_stage = model.get_stage_implicit_inputs(stage)
    if not inputs_for_stage or len(inputs_for_stage) == 0:
        inputs_for_stage = [None]

    modules_in_stage = model.get_modules_by_stage(stage)
    for module_id in modules_in_stage:
        module = modules_in_stage[module_id]
        parameters = model.get_module_parameters(module)
        if not parameters or len(parameters) == 0:
            parameters = [None]

        for param in parameters:
            for inputs in inputs_for_stage:
                latest_stage = None
                if inputs:
                    stages_by_output = list(
                        set(
                            [
                                stage.id
                                for input_id in inputs
                                for stage in [model.get_stage_by_output(input_id)]
                                if stage is not None
                            ]
                        )
                    )
                    latest_stage = sorted(stages_by_output, key=stage_ordering.index)[
                        -1
                    ]
                    explicit_inputs = model.get_explicit_inputs(inputs)
                    inputs = {
                        k: v.replace("{input}", "{pre}")
                        for k, v in explicit_inputs.items()
                    }

                node = BenchmarkNode(
                    model,
                    stage,
                    module,
                    param,
                    inputs,
                    outputs,
                    after=latest_stage,
                )
                nodes.append(node)

    return nodes


def build_benchmark_dag(model: Benchmark, out_dir: Path) -> DiGraph:
    g = DiGraph()

    G_stages = build_stage_dag(model)
    stage_ordering = compute_stage_order(G_stages)

    stage_nodes_map = {}
    for stage_id, stage in model.get_stages().items():
        nodes = expand_stage_nodes(model, stage, out_dir, stage_ordering)
        nodes_with_stage = [(node, {"stage": stage_id}) for node in nodes]
        g.add_nodes_from(nodes_with_stage)
        stage_nodes_map[stage_id] = nodes

    for current_node in g.nodes:
        after_stage = current_node.after
        if after_stage:
            departure_stage_nodes = stage_nodes_map[after_stage]
            for departure_stage_node in departure_stage_nodes:
                g.add_edge(departure_stage_node, current_node)

    return g


def build_stage_dag(model: Benchmark) -> DiGraph:
    g = DiGraph()

    for stage_id, stage in model.get_stages().items():
        g.add_node(stage_id)
        input_ids = [
            input_id for input in (stage.inputs or []) for input_id in input.entries
        ]
        dep_stages = [model.get_output_stage(input_id) for input_id in input_ids]
        for dep in dep_stages:
            if dep is not None:
                g.add_edge(dep.id, stage.id)

    return g


def find_initial_and_terminal_nodes(
    graph: DiGraph,
) -> Tuple[List[BenchmarkNode], List[BenchmarkNode]]:
    initial_nodes = [node for node, in_degree in graph.in_degree() if in_degree == 0]
    terminal_nodes = [
        node for node, out_degree in graph.out_degree() if out_degree == 0
    ]
    return initial_nodes, terminal_nodes


def list_all_paths(graph: DiGraph, source: BenchmarkNode, target: BenchmarkNode):
    all_paths = list(all_simple_paths(graph, source=source, target=target))
    return all_paths


def contains_all(path, modules):
    path_modules = [node.module_id for node in path]
    return all(module in path_modules for module in modules)


def exclude_paths(paths, path_exclusions):
    updated_paths = []
    for path in paths:
        should_exclude = False
        for module, excluded_modules in path_exclusions.items():
            for excluded_module in excluded_modules:
                if contains_all(path, [module, excluded_module]):
                    should_exclude = True

        if not should_exclude:
            updated_paths.append(path)

    return updated_paths


def compute_stage_order(stage_dag: DiGraph) -> List:
    try:
        topological_order = list(topological_sort(stage_dag))

    except NetworkXUnfeasible:
        raise ValidationError(
            "The stage graph has a cyclic dependencies. This benchmark can not be resolved"
        )

    return topological_order


def find_node_by_id(graph, node_id: str) -> Optional[BenchmarkNode]:
    """Find a node in the graph by its ID.

    Args:
        graph: The benchmark DAG graph
        node_id: The ID of the node to find

    Returns:
        The node with matching ID, or None if not found
    """
    for node in graph.nodes:
        if node.get_id() == node_id:
            return node
    return None


def get_nodes_by_module_id(graph, module_id: str) -> List[BenchmarkNode]:
    """Get all nodes in the graph that belong to a specific module.

    Args:
        graph: The benchmark DAG graph
        module_id: The module ID to filter by

    Returns:
        List of nodes belonging to the specified module
    """
    nodes = []
    for node in graph.nodes:
        if node.module_id == module_id:
            nodes.append(node)
    return nodes


def get_nodes_by_stage_id(graph, stage_id: str) -> List[BenchmarkNode]:
    """Get all nodes in the graph that belong to a specific stage.

    Args:
        graph: The benchmark DAG graph
        stage_id: The stage ID to filter by

    Returns:
        List of nodes belonging to the specified stage
    """
    nodes = []
    for node in graph.nodes:
        if node.stage_id == stage_id:
            nodes.append(node)
    return nodes


def find_node_with_module_id(graph, module_id: str) -> Optional[BenchmarkNode]:
    """Find the first node with a specific module ID.

    Args:
        graph: The benchmark DAG graph
        module_id: The module ID to search for

    Returns:
        The first node with matching module ID, or None if not found
    """
    for node in graph.nodes:
        if node.module_id == module_id:
            return node
    return None


def get_benchmark_datasets(model, stages) -> List[str]:
    """Extract dataset IDs from initial stages.

    Args:
        model: The benchmark model
        stages: Dictionary of stage_id -> stage mappings

    Returns:
        List of dataset IDs from initial stages
    """
    datasets = []
    for _, stage in stages.items():
        # There should be only one initial stage
        # TODO(ben): worthy to move this assumption to model validation.
        if model.is_initial(stage):
            for module in stage.modules:
                datasets.append(module.id)
            break
    return datasets
