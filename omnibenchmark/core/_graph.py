from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple, Optional

from omnibenchmark.dag import (
    DiGraph,
    topological_sort,
    all_simple_paths,
    NetworkXUnfeasible,
)
from omnibenchmark.model import Benchmark, Stage, ValidationError

from ._node import BenchmarkNode
from ._paths import is_lineage_excluded


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


def stage_adjacency(model: Benchmark) -> List[Tuple[str, str, List[str]]]:
    """Declared stage-to-stage topology edges.

    Returns one ``(upstream_stage_id, stage_id, shared_output_ids)`` triple per
    pair of stages where the downstream stage's declared inputs reference one or
    more of the upstream stage's outputs. ``shared_output_ids`` is exactly those
    outputs (sorted) — suitable for edge labels. Triples are ordered by
    downstream stage, then by first appearance of the upstream stage in the
    downstream stage's inputs.

    This is the single source of truth for topology edges. It reads the model's
    *declared* inputs directly (every upstream stage of an input group, not just
    the topologically-latest one the execution DAG keeps for output-path
    chaining), so it stays correct when an input group spans multiple stages.

    Consumers:
      * ``upstream_stage_ids`` / ``build_topology_dag`` here, which back the
        mermaid and dot exporters.
      * **obeditor** — the browser-based benchmark editor loads omnibenchmark as
        a Pyodide wheel and calls this to render its pipeline-topology diagram
        (stage-level, labeled edges). It is an *external* consumer, so this
        function has no in-repo caller beyond the wrappers below; do not flag it
        as dead code or inline it away.
    """
    adjacency: List[Tuple[str, str, List[str]]] = []
    for stage in model.get_stages().values():
        shared_by_stage: Dict[str, List[str]] = {}
        order: List[str] = []
        for entries in model.get_stage_implicit_inputs(stage):
            for output_id in entries:
                producing_stage = model.get_output_stage(output_id)
                if producing_stage is None or producing_stage.id == stage.id:
                    continue
                if producing_stage.id not in shared_by_stage:
                    shared_by_stage[producing_stage.id] = []
                    order.append(producing_stage.id)
                if output_id not in shared_by_stage[producing_stage.id]:
                    shared_by_stage[producing_stage.id].append(output_id)
        for up_stage_id in order:
            adjacency.append(
                (up_stage_id, stage.id, sorted(shared_by_stage[up_stage_id]))
            )
    return adjacency


def upstream_stage_ids(model: Benchmark, stage: "Stage") -> List[str]:
    """IDs of every stage feeding *stage* via its declared inputs.

    Thin wrapper over :func:`stage_adjacency`. Unlike the execution DAG's
    ``after`` join (which keeps only the topologically-latest producing stage so
    output-path construction stays a linear ancestry), this returns *all*
    upstream stages — every declared dependency, regardless of how many stages a
    single input group spans.
    """
    return [
        up_stage_id
        for up_stage_id, stage_id, _ in stage_adjacency(model)
        if stage_id == stage.id
    ]


def build_topology_dag(model: Benchmark, graph: DiGraph) -> DiGraph:
    """Build a display-only DAG with every declared cross-stage edge.

    Reuses the nodes of the execution *graph* but recomputes edges from the
    model's declared inputs: every node of each upstream stage connects to every
    node of the consuming stage. The execution DAG connects each node only to
    its single ``after`` stage (correct for linear output-path chaining but
    lossy when an input group spans multiple upstream stages); this restores the
    dropped edges for visualization without perturbing execution.
    """
    display = DiGraph()

    nodes_by_stage: dict = defaultdict(list)
    for node in graph.nodes:
        nodes_by_stage[node.stage_id].append(node)
        display.add_node(node, stage=node.stage_id)

    for stage in model.get_stages().values():
        targets = nodes_by_stage.get(stage.id, [])
        for up_stage_id in upstream_stage_ids(model, stage):
            for source in nodes_by_stage.get(up_stage_id, []):
                for target in targets:
                    display.add_edge(source, target)

    return display


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


def exclude_paths(paths, path_exclusions):
    return [
        path
        for path in paths
        if not is_lineage_excluded({node.module_id for node in path}, path_exclusions)
    ]


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
