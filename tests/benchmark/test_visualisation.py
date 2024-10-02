import re
import random
import networkx as nx
from pathlib import Path

import pytest

from omni.benchmark import Benchmark, dag


def test_plot_computational_graph():
    benchmark_file = "../data/Benchmark_001.yaml"
    benchmark_file_path = Path(__file__).parent / benchmark_file
    output_file = "test_plot_computational_graph.png"
    output_file_path = Path(output_file)

    try:
        with open(benchmark_file_path, "r") as file:
            benchmark = Benchmark(benchmark_file_path)
            benchmark.plot_computational_graph(output_file=output_file)
            assert (
                output_file_path.exists()
            ), f"Output file {output_file} was not created."
    except Exception as e:
        pytest.fail(f"Plotting benchmark computational graph failed: {e}")
    finally:
        if output_file_path.exists():
            output_file_path.unlink()


@pytest.mark.timeout(60)
def test_plot_computational_graph_scaling():
    Glarge = generate_graph(100)
    fig = dag.export_to_figure(
        Glarge,
        title="Large Graph (100 nodes)",
    )
    fig.savefig("large_graph.png")

    Gmedium = generate_graph(50)
    fig = dag.export_to_figure(
        Gmedium,
        title="Medium Graph (50 nodes)",
    )
    fig.savefig("medium_graph.png")

    Gsmall = generate_graph(10)
    fig = dag.export_to_figure(
        Gsmall,
        title="Small Graph (10 nodes)",
    )
    fig.savefig("small_graph.png")


def test_export_to_mermaid():
    benchmark_file = "../data/Benchmark_001.yaml"
    benchmark_file_path = Path(__file__).parent / benchmark_file
    expected_mermaid = """---
    title: Benchmark_001
    ---
    flowchart LR
        \tclassDef param fill:#f96
        \tsubgraph data
            \t\tD1
            \t\tD2
        \tend
        \tsubgraph process
            \t\tP1
            \t\tD1 --> P1
            \t\tD2 --> P1
            \t\tP2
            \t\tD1 --> P2
            \t\tD2 --> P2
        \tend
        \tsubgraph methods
            \t\tM1
            \t\tD1 --> M1
            \t\tD2 --> M1
            \t\tP1 --> M1
            \t\tP2 --> M1
            \t\tM2
            \t\tD1 --> M2
            \t\tD2 --> M2
            \t\tP1 --> M2
            \t\tP2 --> M2
        \tend
        \tsubgraph metrics
            \t\tm1
            \t\tM1 --> m1
            \t\tM2 --> m1
            \t\tm2
            \t\tM1 --> m2
            \t\tM2 --> m2
            \t\tm3
            \t\tM1 --> m3
            \t\tM2 --> m3
        \tend
    """

    try:
        with open(benchmark_file_path, "r") as file:
            benchmark = Benchmark(benchmark_file_path)
            mermaid = benchmark.export_to_mermaid()
            assert clean(mermaid).startswith(clean(expected_mermaid))
    except Exception as e:
        pytest.fail(f"Exporting benchmark topology to mermaid failed: {e}")


def clean(output: str) -> str:
    output = output.strip()
    output = output.replace("    ", "")

    # Replace different newline characters with a single '\n'
    normalized_output = re.sub(r"\r\n|\r", "\n", output)

    # Replace multiple spaces and tabs with a single space
    normalized_output = re.sub(r"[ \t]+", " ", normalized_output)

    return normalized_output


def generate_graph(n_nodes, stage_proportions=None):
    """Generate the graph using the Barabási–Albert model"""

    # Calculate the number of nodes for each stage
    if stage_proportions is None:
        stage_proportions = {
            "datasets": 0.2,  # 20% of the nodes for datasets
            "preprocess": 0.3,  # 30% of the nodes for preprocessing
            "methods": 0.4,  # 40% of the nodes for methods
            "metrics": 0.1,  # 10% of the nodes for metrics
        }

    n_datasets = int(n_nodes * stage_proportions["datasets"])
    n_preprocess = int(n_nodes * stage_proportions["preprocess"])
    n_methods = int(n_nodes * stage_proportions["methods"])
    n_metrics = int(n_nodes * stage_proportions["metrics"])

    # Create a directed graph
    G = nx.DiGraph()

    # Add dataset nodes
    datasets = [f"Dataset_{i + 1}" for i in range(n_datasets)]
    datasets_with_stage = [(d, {"stage": "datasets"}) for d in datasets]
    G.add_nodes_from(datasets_with_stage)

    # Add preprocessing nodes
    preprocesses = [f"Preprocess_{i + 1}" for i in range(n_preprocess)]
    preprocesses_with_stage = [(p, {"stage": "preprocess"}) for p in preprocesses]
    G.add_nodes_from(preprocesses_with_stage)

    # Add methods nodes
    methods = [f"Method_{i + 1}" for i in range(n_methods)]
    methods_with_stage = [(m, {"stage": "methods"}) for m in methods]
    G.add_nodes_from(methods_with_stage)

    # Add metrics nodes
    metrics = [f"Metric_{i + 1}" for i in range(n_metrics)]
    metrics_with_stage = [(m, {"stage": "metrics"}) for m in metrics]
    G.add_nodes_from(metrics_with_stage)

    # Connect datasets to preprocessing steps
    for dataset in datasets:
        for preprocess in preprocesses:
            G.add_edge(dataset, preprocess)

    # Connect preprocessing steps to methods
    for preprocess in preprocesses:
        for method in methods:
            G.add_edge(preprocess, method)

    # Connect methods to metrics
    for method in methods:
        for metric in metrics:
            G.add_edge(method, metric)

    # Generate random connections between stages
    # Randomly connect datasets to methods
    for _ in range(
        n_datasets * n_methods // 2
    ):  # Random connections between datasets and methods
        dataset = random.choice(datasets)
        method = random.choice(methods)
        G.add_edge(dataset, method)

    # Randomly connect preprocess steps to metrics
    for _ in range(
        n_preprocess * n_metrics // 2
    ):  # Random connections between preprocess and metrics
        preprocess = random.choice(preprocesses)
        metric = random.choice(metrics)
        G.add_edge(preprocess, metric)

    return G
