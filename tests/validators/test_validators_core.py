"""Unit tests for the pure helpers in omnibenchmark.core.validators."""

import pytest

from pathlib import Path

from omnibenchmark.core.validators import (
    build_targets,
    discover_validators,
    path_to_glob,
)


@pytest.mark.short
@pytest.mark.parametrize(
    "template,expected",
    [
        ("{dataset}_pcas.tsv", "*_pcas.tsv"),
        ("{dataset}.h5ad", "*.h5ad"),
        ("{dataset}_{params}_out.txt", "*_*_out.txt"),
        ("static.json", "static.json"),
    ],
)
def test_path_to_glob(template, expected):
    assert path_to_glob(template) == expected


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
def test_build_targets_matches_produced_files_to_validators():
    validators = {
        ("PCA", "pcas_tsv"): ("Rscript", "/x/validate.R"),
        ("DATA", "rawdata_h5ad"): ("python3", "/y/validate.py"),
    }
    # get_stage_outputs returns the full expanded path; only the leaf is matched
    stage_outputs = {
        "PCA": {
            "pcas_tsv": "{input}/{stage}/{module}/{params}/{dataset}_pcas.tsv",
            "loadings_tsv": "{input}/{stage}/{module}/{params}/{dataset}_loadings.tsv",
        },
        "DATA": {"rawdata_h5ad": "{input}/{stage}/{module}/{dataset}.h5ad"},
    }
    # produced files are deeply nested under hidden .{param} dirs (real layout)
    produced = [
        ("DATA", "out/DATA/data/.fa42/data.h5ad"),
        ("PCA", "out/DATA/data/.fa42/FILT/f/.default/PCA/p/.77f7/data_pcas.tsv"),
        ("PCA", "out/DATA/data/.fa42/FILT/f/.default/PCA/p/.77f7/data_loadings.tsv"),
    ]
    targets = build_targets(produced, validators, stage_outputs, Path("out"))

    assert targets["DATA/data/.fa42/data.h5ad"] == ["python3", "/y/validate.py"]
    assert targets["DATA/data/.fa42/FILT/f/.default/PCA/p/.77f7/data_pcas.tsv"] == [
        "Rscript",
        "/x/validate.R",
    ]
    # loadings.tsv has no validator (no pcas match) -> not a target
    assert len(targets) == 2
