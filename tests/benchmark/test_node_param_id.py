from pathlib import Path

import pytest

from omnibenchmark.core import BenchmarkExecution


@pytest.mark.short
def test_parameterless_node_uses_dot_default():
    """`ob run` writes parameterless modules under `.default`; status must look there."""
    b = BenchmarkExecution(Path(__file__).parent / "../data/Benchmark_001.yaml")
    node = next(
        n for n in b.get_nodes() if n.stage_id == "data" and n.module_id == "D1"
    )
    assert node.param_id == ".default"
    assert node.get_id() == "data-D1-.default"
