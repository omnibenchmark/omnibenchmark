from pathlib import Path
from omnibenchmark.benchmark import Benchmark

def load(benchmark_yaml: str):
    return Benchmark(Path(benchmark_yaml))


def load_node(benchmark_yaml: str, node_id: str):
    benchmark = load(benchmark_yaml)
    return benchmark.get_node_by_id(node_id)
