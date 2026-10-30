"""Unit tests for the pure helpers in omnibenchmark.core.validators."""

import pytest

from pathlib import Path

from omnibenchmark.core.validators import (
    build_targets,
    discover_validators,
)


@pytest.mark.short
def test_discover_validators(tmp_path):
    # keyed by (stage_id, output_id), e.g. PCA/pcas_tsv
    d = tmp_path / "validators" / "PCA" / "pcas_tsv"
    d.mkdir(parents=True)
    (d / "validate.R").write_text("# noop")
    found = discover_validators(tmp_path / "validators")
    assert set(found) == {("PCA", "pcas_tsv")}
    interp, script = found[("PCA", "pcas_tsv")]
    assert interp == "Rscript"
    assert script.endswith("validate.R")


@pytest.mark.short
def test_discover_validators_rejects_multiple(tmp_path):
    d = tmp_path / "validators" / "s" / "o"
    d.mkdir(parents=True)
    (d / "validate.py").write_text("")
    (d / "validate.sh").write_text("")
    with pytest.raises(ValueError):
        discover_validators(tmp_path / "validators")


@pytest.mark.short
def test_build_targets_keys_by_output_id():
    validators = {
        ("PCA", "pcas_tsv"): ("Rscript", "/x/validate.R"),
        ("DATA", "rawdata_h5ad"): ("python3", "/y/validate.py"),
    }
    produced = [
        ("DATA", "rawdata_h5ad", "out/data.data/data.h5ad"),
        ("PCA", "pcas_tsv", "out/data.data/pca.p/data_pcas.tsv"),
        # no validator for loadings -> not a target
        ("PCA", "loadings_tsv", "out/data.data/pca.p/data_loadings.tsv"),
    ]
    assert build_targets(produced, validators, Path("out")) == {
        "data.data/data.h5ad": ["python3", "/y/validate.py"],
        "data.data/pca.p/data_pcas.tsv": ["Rscript", "/x/validate.R"],
    }


@pytest.mark.short
def test_validators_require_api_0_8():
    from omnibenchmark.model.benchmark import Benchmark

    plan = """\
id: b
version: "1.0"
benchmarker: t
api_version: "0.7.0"
software_backend: host
software_environments:
  host:
    description: host
stages:
  - id: data
    modules:
      - id: D1
        software_environment: host
        repository: {url: https://example.org/x.git, commit: abc1234}
    outputs:
      - id: raw
        path: "{module.id}_data.json"
validators:
  env: host
"""
    with pytest.raises(Exception, match="requires api_version"):
        Benchmark.from_yaml(plan)
