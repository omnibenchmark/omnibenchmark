from pathlib import Path
from typing import List, Tuple

import networkx as nx
import pydot

import omni_schema.datamodel.omni_schema

from omnibenchmark.benchmark.converter import LinkMLConverter
from omnibenchmark.benchmark.benchmark_node import BenchmarkNode
from omnibenchmark.benchmark.validation import ValidationError


def expend_stage_nodes(
    converter: LinkMLConverter,
    stage: omni_schema.datamodel.omni_schema.Stage,
    out_dir: Path,
    stage_ordering: List[str],
) -> List[BenchmarkNode]:
    nodes = []

    input_dirname = str(out_dir) if converter.is_initial(stage) else "{pre}"
    stage_outputs = converter.get_stage_outputs(stage).values()
    outputs = [x.replace("{input}", input_dirname) for x in stage_outputs]

    inputs_for_stage = converter.get_stage_implicit_inputs(stage)
    if not inputs_for_stage or len(inputs_for_stage) == 0:
        inputs_for_stage = [None]

    modules_in_stage = converter.get_modules_by_stage(stage)
    for module_id in modules_in_stage:
        module = modules_in_stage[module_id]
        parameters = converter.get_module_parameters(module)
        if not parameters or len(parameters) == 0:
            parameters = [None]

        for param in parameters:
            for inputs in inputs_for_stage:
                latest_stage = None
                if inputs:
                    stages_by_output = list(
                        set(
                            [
                                converter.get_stage_by_output(input_id).id
                                for input_id in inputs
                            ]
                        )
                    )
                    latest_stage = sorted(stages_by_output, key=stage_ordering.index)[
                        -1
                    ]
                    explicit_inputs = converter.get_explicit_inputs(inputs)
                    inputs = {
                        k: v.replace("{input}", "{pre}")
                        for k, v in explicit_inputs.items()
                    }

                node = BenchmarkNode(
                    converter,
                    stage,
                    module,
                    param,
                    inputs,
                    outputs,
                    after=latest_stage,
                )
                nodes.append(node)

    return nodes


def build_benchmark_dag(converter: LinkMLConverter, out_dir: Path) -> nx.DiGraph:
    g = nx.DiGraph()

    G_stages = build_stage_dag(converter)
    stage_ordering = compute_stage_order(G_stages)

    stage_nodes_map = {}
    for stage_id, stage in converter.get_stages().items():
        nodes = expend_stage_nodes(converter, stage, out_dir, stage_ordering)
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


def build_stage_dag(converter: LinkMLConverter) -> nx.DiGraph:
    g = nx.DiGraph()

    for stage_id, stage in converter.get_stages().items():
        g.add_node(stage_id)
        input_ids = [input_id for input in stage.inputs for input_id in input.entries]
        dep_stages = [converter.get_output_stage(input_id) for input_id in input_ids]
        for dep in dep_stages:
            g.add_edge(dep.id, stage.id)

    return g


def find_initial_and_terminal_nodes(
    graph: nx.DiGraph,
) -> Tuple[List[BenchmarkNode], List[BenchmarkNode]]:
    initial_nodes = [node for node, in_degree in graph.in_degree() if in_degree == 0]
    terminal_nodes = [
        node for node, out_degree in graph.out_degree() if out_degree == 0
    ]
    return initial_nodes, terminal_nodes


def list_all_paths(graph: nx.DiGraph, source: BenchmarkNode, target: BenchmarkNode):
    all_paths = list(nx.all_simple_paths(graph, source=source, target=target))
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


def compute_stage_order(stage_dag: nx.DiGraph) -> List:
    try:
        topological_order = list(nx.topological_sort(stage_dag))

    except nx.NetworkXUnfeasible:
        raise ValidationError(
            "The stage graph has a cyclic dependencies. This benchmark can not be resolved"
        )

    return topological_order


def export_to_dot(
    G: nx.DiGraph,
    title: str = None,
):
    import matplotlib.pyplot as plt

    # Dynamically scale the node size based on node count
    nodes_count = len(G.nodes)
    div_nodes_count = max(1, nodes_count // 10)
    graph_size = max(15, 15 * div_nodes_count)

    # Color nodes by stage (assuming 'stage' is a node attribute)
    stages = nx.get_node_attributes(G, "stage", default="none")
    unique_stages = list(set(stages.values()))  # Get unique stages

    # Define a colormap with different shades for the stages
    stage_colors = plt.get_cmap("inferno", max(len(unique_stages), 5))

    # Convert the graph to a PyDot graph object
    pydot_graph = pydot.Dot(
        graph_type="digraph", strict=True, label=title, labelloc="top", fontsize=20
    )
    pydot_graph.set_graph_defaults(
        size=f"{graph_size},{graph_size}!", ratio="fill", margin=div_nodes_count
    )

    # Define the style for nodes
    node_defaults = {
        "shape": "rect",
        "style": "filled,rounded",
        "fontsize": "12",
        "fontcolor": "white",
        "width": "1.5",
        "height": "0.6",
        "penwidth": "1.0",
    }
    pydot_graph.set_node_defaults(**node_defaults)

    # Define the style for edges
    edge_defaults = {"color": "#CCCCCC", "penwidth": "0.5", "arrowsize": "0.7"}
    pydot_graph.set_edge_defaults(**edge_defaults)

    for node in G.nodes:
        node_name = str(node)
        rgba_color = stage_colors(unique_stages.index(stages[node]))
        hex_color = _rgba_to_hex(rgba_color)
        pydot_node = pydot.Node(node_name, label=node_name, fillcolor=hex_color)
        pydot_graph.add_node(pydot_node)

    for source, target in G.edges:
        pydot_edge = pydot.Edge(str(source), str(target))
        pydot_graph.add_edge(pydot_edge)

    return pydot_graph


def _rgba_to_hex(rgba):
    r = int(rgba[0] * 255)
    g = int(rgba[1] * 255)
    b = int(rgba[2] * 255)
    return f"#{r:02X}{g:02X}{b:02X}"
