from pathlib import Path
from omnibenchmark.benchmark import Benchmark

def load(benchmark_yaml: str, config):
    out_dir = config['out_dir']
    return Benchmark(Path(benchmark_yaml), Path(out_dir))


def load_node(benchmark_yaml: str, node_id: str):
    benchmark = Benchmark(Path(benchmark_yaml), Path("/tmp"))
    return benchmark.get_node_by_id(node_id)
