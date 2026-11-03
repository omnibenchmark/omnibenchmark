"""Unit tests for the pure helpers in omnibenchmark.core.validators."""

import pytest

from pathlib import Path

from omnibenchmark.core.validators import (
    build_targets,
    collect_validators,
    discover_validators,
    verdict_status,
)
from omnibenchmark.model.benchmark import Benchmark


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
    # missing scripts hash to "" (the rule then records exit 127)
    assert build_targets(produced, validators, Path("out")) == {
        "data.data/data.h5ad": ["python3", "/y/validate.py", ""],
        "data.data/pca.p/data_pcas.tsv": ["Rscript", "/x/validate.R", ""],
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


PLAN = """\
id: b
version: "1.0"
benchmarker: t
api_version: "0.8.0"
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
      - id: meta
        path: "{module.id}_meta.json"
        validator:
          script: checks/meta.sh
validators: {}
"""


@pytest.mark.short
def test_collect_validators_declared_wins_and_typos_reported(tmp_path):
    vdir = tmp_path / "validators"
    for key in ("data/raw", "data/meta", "data/rwa"):
        (vdir / key).mkdir(parents=True)
        (vdir / key / "validate.py").write_text("")
    found, unmatched = collect_validators(Benchmark.from_yaml(PLAN), tmp_path, vdir)
    assert found[("data", "raw")][0] == "python3"
    assert found[("data", "meta")] == ("bash", str(tmp_path / "checks" / "meta.sh"))
    assert unmatched == {("data", "rwa")}


@pytest.mark.short
def test_build_targets_hashes_script(tmp_path):
    script = tmp_path / "v.py"
    script.write_text("pass\n")
    targets = build_targets(
        [("s", "o", "out/a.txt")], {("s", "o"): ("python3", str(script))}, Path("out")
    )
    sha = targets["a.txt"][2]
    script.write_text("raise SystemExit(1)\n")
    targets = build_targets(
        [("s", "o", "out/a.txt")], {("s", "o"): ("python3", str(script))}, Path("out")
    )
    assert len(sha) == 64 and targets["a.txt"][2] != sha


@pytest.mark.short
def test_verdict_status():
    assert verdict_status(0) == "pass"
    assert verdict_status(1) == "fail"
    assert verdict_status(127) == "error"
    assert verdict_status(126) == "error"
    assert verdict_status(None) == "error"


@pytest.mark.short
def test_validators_env_must_be_declared():
    with pytest.raises(Exception, match="for validators is not declared"):
        Benchmark.from_yaml(PLAN.replace("validators: {}", "validators: {env: nope}"))


@pytest.mark.short
def test_declared_validator_script_must_exist(tmp_path):
    model = Benchmark.from_yaml(PLAN)
    with pytest.raises(Exception, match="checks/meta.sh"):
        model.validate_execution_context(tmp_path)
    (tmp_path / "checks").mkdir()
    (tmp_path / "checks" / "meta.sh").write_text("")
    model.validate_execution_context(tmp_path)


@pytest.mark.short
def test_declared_validator_requires_api_0_8():
    plan = PLAN.replace('"0.8.0"', '"0.7.0"').replace("validators: {}\n", "")
    with pytest.raises(Exception, match="requires api_version"):
        Benchmark.from_yaml(plan)
