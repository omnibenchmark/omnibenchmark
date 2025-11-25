from omnibenchmark.dag import DiGraph, get_node_attributes


def export_to_dot(
    G: DiGraph,
    title: str = "",
):
    import matplotlib.pyplot as plt
    import pydot

    # Dynamically scale the node size based on node count
    nodes_count = len(G.nodes)
    div_nodes_count = max(1, nodes_count // 10)
    graph_size = max(15, 15 * div_nodes_count)

    # Color nodes by stage (assuming 'stage' is a node attribute)
    stages = get_node_attributes(G, "stage", default="none")
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
