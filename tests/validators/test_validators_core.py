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
