"""Output paths a node predicts for `ob describe status` / `ob validate outputs`."""

import pytest

from omnibenchmark.core import BenchmarkExecution

_PLAN = """
id: t
description: t
version: "1.0"
benchmarker: t
api_version: "{api}"
software_backend: host
software_environments:
  host:
    description: host
stages:
  - id: data
    modules:
      - id: D1
        software_environment: host
        repository:
          url: https://example.org/data.git
          commit: abc1234
    outputs:
      - id: counts
        path: "counts.txt"
"""


@pytest.mark.short
@pytest.mark.parametrize(
    "api, expected",
    [
        # 0.8.0 writes no `.default` for a node without parameters
        ("0.8.0", "out/data.D1/counts.txt"),
        ("0.7.0", "out/data/D1/.default/counts.txt"),
    ],
)
def test_parameterless_node_path(tmp_path, api, expected):
    plan = tmp_path / "plan.yaml"
    plan.write_text(_PLAN.format(api=api))
    benchmark = BenchmarkExecution(plan, tmp_path / "out")
    (node,) = benchmark.get_nodes()
    paths = node.get_output_paths({"input": "out", "dataset": "D1"})
    (path,) = paths
    assert path == str(tmp_path / expected)
