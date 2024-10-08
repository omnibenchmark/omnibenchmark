import math
import random
from typing import List, Tuple

import networkx as nx
import matplotlib.pyplot as plt
import omni_schema.datamodel.omni_schema
from matplotlib.figure import Figure

from omni.benchmark.converter import LinkMLConverter
from omni.benchmark.benchmark_node import BenchmarkNode
from omni.benchmark.validation import ValidationError
from omni.constants import LayoutDesign


def expend_stage_nodes(
    converter: LinkMLConverter,
    stage: omni_schema.datamodel.omni_schema.Stage,
    out_dir: str,
    stage_ordering: List[str],
) -> List[BenchmarkNode]:
    nodes = []

    input_dirname = out_dir if converter.is_initial(stage) else "{pre}"
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

        for param_id, param in enumerate(parameters):
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
                    param_id,
                    after=latest_stage,
                )
                nodes.append(node)

    return nodes


def build_benchmark_dag(converter: LinkMLConverter, out_dir: str) -> nx.DiGraph:
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


def export_to_figure(
    G: nx.DiGraph,
    layout_design: LayoutDesign = LayoutDesign.Hierarchical,
    title: str = None,
) -> Figure:
    if layout_design == LayoutDesign.Spring:
        layout = nx.circular_layout(G)
    elif layout_design == LayoutDesign.Hierarchical:
        from networkx.drawing.nx_agraph import pygraphviz_layout

        layout = pygraphviz_layout(G, prog="neato")
    else:
        raise ValueError(
            f"Graph can only be exported using the ${LayoutDesign.value} layouts."
        )

    # Dynamically scale the figure size based on node count
    nodes_count = len(G.nodes)
    figure_size = (max(15, nodes_count // 4), max(15, nodes_count // 4))

    # Create a new figure and set the size
    fig = plt.figure(figsize=figure_size)
    plt.gca().set_facecolor("white")

    # Draw edges with arrows to show flow direction
    arrow_size = 20 / math.log10(nodes_count)
    nx.draw_networkx_edges(
        G,
        layout,
        arrowstyle="-|>",
        arrowsize=arrow_size,
        edge_color="#AAAAAA",
        alpha=0.6,
        width=0.5,
    )

    # Color nodes by stage (assuming 'stage' is a node attribute)
    stages = nx.get_node_attributes(
        G, "stage", default="none"
    )  # Node attribute 'stage'
    unique_stages = list(set(stages.values()))  # Get unique stages

    # Define a colormap with different shades for the stages
    stage_colors = plt.get_cmap(
        "inferno", max(len(unique_stages), 5)
    )  # Use Set3 colormap for distinct stage colors
    node_colors = [stage_colors(unique_stages.index(stages[n])) for n in G.nodes()]

    node_size = 500 / math.log10(nodes_count)
    node_sizes = [node_size] * nodes_count
    nx.draw_networkx_nodes(
        G,
        layout,
        nodelist=G.nodes(),
        node_size=node_sizes,
        node_color=node_colors,
    )

    # Create adjusted label positions
    # Vertical offset to position below the node
    labels_pos = {}

    for n in layout:
        # Add some random vertical variation
        if layout_design == LayoutDesign.Hierarchical:
            label_offset_y = -node_size / 100
            y_offset = label_offset_y - random.uniform(0, 5)
        else:
            label_offset_y = -node_size / 5000
            y_offset = label_offset_y

        labels_pos[n] = (layout[n][0], layout[n][1] + y_offset)

    # Draw the labels with adjusted positions
    nx.draw_networkx_labels(
        G,
        labels_pos,
        labels={n: n for n in G.nodes()},
        font_size=8,
        font_color="black",
        bbox=dict(facecolor="white", edgecolor="black", boxstyle="round,pad=0.2"),
    )

    # Set the title if provided
    if title:
        plt.title(title, fontsize=14)

    return fig
