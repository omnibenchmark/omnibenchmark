from pathlib import Path

from omni.benchmark import Benchmark


def test_parse_benchmark():
    benchmark_file = "../data/Benchmark_001.yaml"
    benchmark_file_path = Path(__file__).parent / benchmark_file

    benchmark = Benchmark(benchmark_file_path)
