import re
import random

import pytest
from pathlib import Path

from omnibenchmark.core import Benchmark
from omnibenchmark.core._dot import export_to_dot
from omnibenchmark.core._graph import (
    build_benchmark_dag,
    build_topology_dag,
    stage_adjacency,
)
from omnibenchmark.model import Benchmark as BenchmarkModel
from omnibenchmark.dag import DiGraph


def _module_edges(model):
    """Module-level edges of the display (topology) DAG built from *model*."""
    exec_dag = build_benchmark_dag(model, Path("."))
    display = build_topology_dag(model, exec_dag)
    return {(u.module_id, v.module_id) for u, v in display.edges}


_TWO_STAGE_DIAMOND = """\
id: diamond
description: d
version: "1.0"
benchmarker: x
storage: https://s/
api_version: "0.3.0"
storage_api: S3
storage_bucket_name: b
software_backend: host
software_environments:
  e: {description: e, easyconfig: e.eb, envmodule: e, conda: e.yaml, apptainer: http://x/y.sif}
stages:
  - id: data
    modules:
      - {id: D, software_environment: e, repository: {url: https://github.com/o/data.git, commit: 63b7b36}}
    outputs:
      - {id: data.out, path: "{dataset}.txt"}
  - id: proc
    modules:
      - {id: P, software_environment: e, repository: {url: https://github.com/o/p.git, commit: 63b7b36}}
    inputs: [data.out]
    outputs:
      - {id: proc.out, path: "p.txt"}
  - id: combine
    modules:
      - {id: C, software_environment: e, repository: {url: https://github.com/o/c.git, commit: 63b7b36}}
    inputs: [data.out, proc.out]
    outputs:
      - {id: combine.out, path: "c.txt"}
"""


def test_topology_dag_restores_cross_stage_edges():
    """The display DAG draws an edge from *every* declared upstream stage, even
    when one input group spans several (regression: only the latest survived)."""
    model = BenchmarkModel.from_yaml(_TWO_STAGE_DIAMOND)
    assert _module_edges(model) == {("D", "P"), ("D", "C"), ("P", "C")}


def test_topology_dag_does_not_mutate_execution_dag():
    """Building the display DAG must not perturb the execution DAG's edges
    (which must stay a linear ancestry for correct output-path construction)."""
    model = BenchmarkModel.from_yaml(_TWO_STAGE_DIAMOND)
    exec_dag = build_benchmark_dag(model, Path("."))
    before = {(u.module_id, v.module_id) for u, v in exec_dag.edges}
    build_topology_dag(model, exec_dag)
    after = {(u.module_id, v.module_id) for u, v in exec_dag.edges}
    assert before == after == {("D", "P"), ("P", "C")}


def test_export_to_dot_contains_cross_stage_edge():
    """The rendered dot graph includes the restored data->combine edge."""
    model = BenchmarkModel.from_yaml(_TWO_STAGE_DIAMOND)
    display = build_topology_dag(model, build_benchmark_dag(model, Path(".")))
    dot = export_to_dot(display, title="diamond")
    rendered = str(dot)
    assert "combine-C" in rendered and "data-D" in rendered


def test_stage_adjacency_edges_and_labels():
    """stage_adjacency yields every declared edge with its shared output ids."""
    model = BenchmarkModel.from_yaml(_TWO_STAGE_DIAMOND)
    adjacency = stage_adjacency(model)
    assert adjacency == [
        ("data", "proc", ["data.out"]),
        ("data", "combine", ["data.out"]),
        ("proc", "combine", ["proc.out"]),
    ]


def test_stage_adjacency_groups_multiple_shared_outputs():
    """Multiple outputs shared between the same stage pair collapse to one edge
    whose label lists all of them (sorted)."""
    model = BenchmarkModel.from_yaml(
        Path(__file__).parent / "../data/Benchmark_001.yaml"
    )
    adjacency = stage_adjacency(model)
    # metrics declares methods.mapping + data.meta + data.data_specific_params
    metrics_from_data = [
        shared for up, down, shared in adjacency if down == "metrics" and up == "data"
    ]
    assert metrics_from_data == [["data.data_specific_params", "data.meta"]]
    # ... and a single edge from methods as well
    assert any(up == "methods" and down == "metrics" for up, down, _ in adjacency)


def test_upstream_stage_ids_consistent_with_adjacency():
    """The wrapper returns exactly the upstream stages from the adjacency."""
    model = BenchmarkModel.from_yaml(_TWO_STAGE_DIAMOND)
    combine = model.get_stage("combine")
    from omnibenchmark.core._graph import upstream_stage_ids

    assert upstream_stage_ids(model, combine) == ["data", "proc"]


@pytest.mark.short
def test_export_computational_to_dot(tmp_path):
    import shutil

    benchmark_file = "../data/Benchmark_001.yaml"
    benchmark_file_path = Path(__file__).parent / benchmark_file
    output_file_path = tmp_path / "computational_graph.dot"

    # Copy benchmark file to tmp_path to avoid writing to current directory
    copied_benchmark_path = tmp_path / "Benchmark_001.yaml"
    shutil.copy(benchmark_file_path, copied_benchmark_path)

    try:
        benchmark = Benchmark(copied_benchmark_path, out_dir=tmp_path / "out")
        dot = benchmark.export_to_dot()
        dot.write(str(output_file_path))
        assert (
            output_file_path.exists()
        ), f"Output file {output_file_path} was not created."
    except Exception as e:
        pytest.fail(f"Plotting benchmark computational graph failed: {e}")


@pytest.mark.timeout(60)
def test_export_computational_to_dot_scaling(tmp_path):
    Glarge = generate_graph(100)
    dot = export_to_dot(
        Glarge,
        title="Large Graph (100 nodes)",
    )
    dot.write(str(tmp_path / "large_graph.dot"))

    Gmedium = generate_graph(50)
    dot = export_to_dot(
        Gmedium,
        title="Medium Graph (50 nodes)",
    )
    dot.write(str(tmp_path / "medium_graph.dot"))

    Gsmall = generate_graph(10)
    dot = export_to_dot(
        Gsmall,
        title="Small Graph (10 nodes)",
    )
    dot.write(str(tmp_path / "small_graph.dot"))


@pytest.mark.short
def test_export_topology_to_mermaid(tmp_path):
    import shutil

    benchmark_file = "../data/Benchmark_001.yaml"
    benchmark_file_path = Path(__file__).parent / benchmark_file

    # Copy benchmark file to tmp_path to avoid writing to current directory
    copied_benchmark_path = tmp_path / "Benchmark_001.yaml"
    shutil.copy(benchmark_file_path, copied_benchmark_path)
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
            \t\tD1 --> m1
            \t\tD2 --> m1
            \t\tM1 --> m1
            \t\tM2 --> m1
            \t\tm2
            \t\tD1 --> m2
            \t\tD2 --> m2
            \t\tM1 --> m2
            \t\tM2 --> m2
            \t\tm3
            \t\tD1 --> m3
            \t\tD2 --> m3
            \t\tM1 --> m3
            \t\tM2 --> m3
        \tend
    """

    try:
        benchmark = Benchmark(copied_benchmark_path, out_dir=tmp_path / "out")
        mermaid = benchmark.export_to_mermaid()
        assert clean(mermaid).startswith(clean(expected_mermaid))
    except Exception as e:
        pytest.fail(f"Exporting benchmark topology to mermaid failed: {e}")


# TODO: get from somewhere else
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
    G = DiGraph()

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
