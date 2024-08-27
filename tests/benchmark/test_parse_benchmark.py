from pathlib import Path

import pytest
from omni_schema.datamodel import omni_schema

from omni.utils import parse_instance


def test_parse_benchmark():
    benchmark_file = "../data/Benchmark_001.yaml"
    benchmark_file_path = Path(__file__).parent / benchmark_file

    try:
        benchmark = parse_instance(benchmark_file_path, omni_schema.Benchmark)
    except Exception as e:
        pytest.fail(f"Parsing benchmark model failed: {e}")
