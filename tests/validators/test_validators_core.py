"""Unit tests for the pure helpers in omnibenchmark.core.validators."""

import pytest

from pathlib import Path

from omnibenchmark.core.validators import (
    build_targets,
    declared_validators,
    missing_scripts,
    stage_scripts,
    verdict_status,
)
from omnibenchmark.model.benchmark import Benchmark, Validator


@pytest.mark.short
def test_build_targets_one_per_file_and_validator():
    validators = {
        ("PCA", "pcas_tsv"): [("Rscript", "/x/a.R"), ("python3", "/x/b.py")],
        ("DATA", "rawdata_h5ad"): [("python3", "/y/validate.py")],
    }
    produced = [
        ("DATA", "rawdata_h5ad", "out/data.data/data.h5ad"),
        ("PCA", "pcas_tsv", "out/data.data/pca.p/data_pcas.tsv"),
        # no validator for loadings -> not a target
        ("PCA", "loadings_tsv", "out/data.data/pca.p/data_loadings.tsv"),
    ]
    pcas = "data.data/pca.p/data_pcas.tsv"
    # missing scripts hash to "" (the rule then records exit 127)
    assert build_targets(produced, validators, Path("out")) == {
        "data.data/data.h5ad.0": [
            "data.data/data.h5ad",
            "python3",
            "/y/validate.py",
            "",
        ],
        f"{pcas}.0": [pcas, "Rscript", "/x/a.R", ""],
        f"{pcas}.1": [pcas, "python3", "/x/b.py", ""],
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
        validators:
          - script: checks/meta.sh
          - script: checks/meta.py
validators: {}
"""


@pytest.mark.short
def test_declared_validators(tmp_path):
    found = declared_validators(Benchmark.from_yaml(PLAN), tmp_path)
    assert found == {
        ("data", "meta"): [
            ("bash", str(tmp_path / "checks" / "meta.sh")),
            ("python3", str(tmp_path / "checks" / "meta.py")),
        ]
    }


@pytest.mark.short
def test_build_targets_hashes_script(tmp_path):
    script = tmp_path / "v.py"
    script.write_text("pass\n")
    targets = build_targets(
        [("s", "o", "out/a.txt")], {("s", "o"): [("python3", str(script))]}, Path("out")
    )
    sha = targets["a.txt.0"][3]
    script.write_text("raise SystemExit(1)\n")
    targets = build_targets(
        [("s", "o", "out/a.txt")], {("s", "o"): [("python3", str(script))]}, Path("out")
    )
    assert len(sha) == 64 and targets["a.txt.0"][3] != sha


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
def test_missing_validator_scripts(tmp_path):
    model = Benchmark.from_yaml(PLAN)
    # A missing optional check must not block `ob run`.
    model.validate_execution_context(tmp_path)
    (tmp_path / "checks").mkdir()
    (tmp_path / "checks" / "meta.sh").write_text("")
    assert missing_scripts(declared_validators(model, tmp_path)) == [
        str(tmp_path / "checks" / "meta.py")
    ]


@pytest.mark.short
def test_metric_collectors_rejected_on_0_8():
    plan = (
        PLAN
        + """\
metric_collectors:
  - id: c
    software_environment: host
    repository: {url: https://example.org/x.git, commit: abc1234}
    inputs: [raw]
    outputs:
      - id: table
        path: "{input}/table.tsv"
"""
    )
    with pytest.raises(Exception, match="removed in api_version 0.8.0"):
        Benchmark.from_yaml(plan)


@pytest.mark.short
def test_declared_validator_requires_api_0_8():
    plan = PLAN.replace('"0.8.0"', '"0.7.0"').replace("validators: {}\n", "")
    with pytest.raises(Exception, match="requires api_version"):
        Benchmark.from_yaml(plan)


# summary_hash() of this plan on main, before validators existed. A model field
# that leaks into the hash changes it and breaks obfilter parents and 003 tags.
TRIMMED = Path(__file__).parents[1] / "data" / "Benchmark_001_trimmed.yaml"
TRIMMED_HASH = "c737a271903debabd86745d7835eed17d8fdac77cb8d5a44ef0c66681e86ef17"


@pytest.mark.short
def test_summary_hash_is_pinned():
    assert Benchmark.from_yaml(TRIMMED.read_text()).summary_hash() == TRIMMED_HASH


@pytest.mark.short
def test_validators_do_not_change_summary_hash():
    b = Benchmark.from_yaml(TRIMMED.read_text())
    b.stages[0].outputs[0].validators = [Validator(script="checks/a.py")]
    assert b.summary_hash() == TRIMMED_HASH


@pytest.mark.short
def test_stage_scripts_copies_by_sha_and_prunes(tmp_path):
    script = tmp_path / "check.py"
    script.write_text("pass\n")
    staged = tmp_path / "out" / ".validation" / "scripts"
    staged.mkdir(parents=True)
    (staged / "old").write_text("gone")
    stage_scripts({"f.0": ["f", "python3", str(script), "abc"]}, tmp_path / "out")
    assert [p.name for p in staged.iterdir()] == ["abc"]
    assert (staged / "abc").read_text() == "pass\n"
