"""Unit tests for exclude-reachability validation (omnibenchmark#289).

`detect_unsatisfiable_excludes()` flags modules whose `exclude` rules leave no
valid upstream lineage, so they would expand to zero nodes at run time.
"""

import textwrap
from pathlib import Path

from omnibenchmark.model import Benchmark

_HEADER = """
id: t
version: "1.0"
benchmarker: "t"
api_version: "0.3.0"
software_backend: host
software_environments:
  py:
    description: py
    envmodule: none
stages:
"""


def _model(tmp_path: Path, stages_yaml: str) -> Benchmark:
    yaml_file = tmp_path / "b.yaml"
    yaml_file.write_text(_HEADER + textwrap.dedent(stages_yaml))
    return Benchmark.from_yaml(yaml_file)


# Two DATA collections; a downstream module per collection, each excluding the
# other collection. Every module still has a valid lineage — nothing is zero-node.
_BENIGN = """
  - id: DATA
    modules:
      - {id: dA, software_environment: py, repository: {url: u, commit: c}}
      - {id: dB, software_environment: py, repository: {url: u, commit: c}}
    outputs:
      - {id: raw, path: "{dataset}.txt"}
  - id: PROC
    inputs: [raw]
    modules:
      - {id: pA, software_environment: py, exclude: [dB], repository: {url: u, commit: c}}
      - {id: pB, software_environment: py, exclude: [dA], repository: {url: u, commit: c}}
    outputs:
      - {id: proc.out, path: "{dataset}.out"}
"""

# FEAT excludes dB (runs only on dA); PCA excludes dA (runs only on dB) but consumes
# FEAT's output, which never exists on dB → PCA is zero-node.
_ZERO_NODE = """
  - id: DATA
    modules:
      - {id: dA, software_environment: py, repository: {url: u, commit: c}}
      - {id: dB, software_environment: py, repository: {url: u, commit: c}}
    outputs:
      - {id: raw, path: "{dataset}.txt"}
  - id: FEAT
    inputs: [raw]
    modules:
      - {id: fe, software_environment: py, exclude: [dB], repository: {url: u, commit: c}}
    outputs:
      - {id: feat.out, path: "{dataset}.feat"}
  - id: PCA
    inputs: [feat.out]
    modules:
      - {id: pc, software_environment: py, exclude: [dA], repository: {url: u, commit: c}}
    outputs:
      - {id: pca.out, path: "{dataset}.pca"}
"""


def test_benign_excludes_are_not_flagged(tmp_path):
    model = _model(tmp_path, _BENIGN)
    assert model.detect_unsatisfiable_excludes() == []


def test_zero_node_module_is_flagged(tmp_path):
    model = _model(tmp_path, _ZERO_NODE)
    errors = model.detect_unsatisfiable_excludes()
    # Only the primary cause (PCA/pc) is reported, not FEAT (which can run).
    assert len(errors) == 1
    assert "PCA/pc" in errors[0]
    assert "feat.out" in errors[0]
    assert "can never run" in errors[0]
