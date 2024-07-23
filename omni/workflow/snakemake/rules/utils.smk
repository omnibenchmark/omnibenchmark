from pathlib import Path
from omni.benchmark import Benchmark

def load(benchmark_yaml: Path):
    return Benchmark(benchmark_yaml)


def load_node(benchmark_yaml: Path, node_id: str):
    benchmark = load(benchmark_yaml)
    return benchmark.get_node_by_id(node_id)