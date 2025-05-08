import os
from pathlib import Path
from omnibenchmark.benchmark import Benchmark

# We are propagating contextual information about the out directory
# via the environment variable that snakemake execution inherits.
# We might want to revisit this to use the snakemake parameters object,
# which might prove to be cleaner in the long run.
OUT_DIR = os.environ.get('OB_OUT_DIR', "out")

def load(benchmark_yaml: str):
    return Benchmark(Path(benchmark_yaml), Path(OUT_DIR))


def load_node(benchmark_yaml: str, node_id: str):
    benchmark = load(benchmark_yaml)
    return benchmark.get_node_by_id(node_id)
