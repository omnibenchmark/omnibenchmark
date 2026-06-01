"""Tests for the `collect performance` CLI command and its helpers."""

import csv
import json
from pathlib import Path

import pytest

from omnibenchmark.cli.collect import (
    _build_record,
    _find_performance_files,
    _read_performance_row,
    _triples,
)

_PERF_HEADER = (
    "s\th:m:s\tmax_rss\tmax_vms\tmax_uss\tmax_pss\tio_in\tio_out\tmean_load\tcpu_time\n"
)
_PERF_ROW = "290.07\t0:04:50\t8031.62\t13551.20\t7742.45\t7812.24\t87.38\t3399.00\t84.93\t247.24\n"


def _write_perf(path: Path, name: str = "performance.txt") -> Path:
    path.mkdir(parents=True, exist_ok=True)
    perf = path / name
    perf.write_text(_PERF_HEADER + _PERF_ROW)
    return perf


@pytest.mark.short
def test_read_performance_row_drops_hms_and_coerces():
    """The h:m:s column is dropped and numbers are coerced to float."""
    import tempfile

    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
        f.write(_PERF_HEADER + _PERF_ROW)
        tmp = Path(f.name)
    try:
        row = _read_performance_row(tmp)
        assert "h:m:s" not in row
        assert row["s"] == 290.07
        assert row["max_rss"] == 8031.62
    finally:
        tmp.unlink()


@pytest.mark.short
def test_triples_groups_path_components():
    parts = ["download", "norman", ".68ee97fd", "split", "simulation", ".a4bbadc1"]
    assert _triples(parts) == [
        ("download", "norman", ".68ee97fd"),
        ("split", "simulation", ".a4bbadc1"),
    ]


@pytest.mark.short
def test_find_performance_files_matches_both_names_and_skips_noise(tmp_path):
    """Bare and dataset-prefixed performance files are found; noise dirs pruned."""
    _write_perf(tmp_path / "download" / "norman" / ".abc", "performance.txt")
    _write_perf(tmp_path / "download" / "old_ds" / ".def", "old_ds_performance.txt")
    # Noise that must be skipped.
    _write_perf(tmp_path / ".modules" / "vendor" / "x", "performance.txt")
    _write_perf(tmp_path / ".snakemake" / "y" / "z", "performance.txt")

    found = _find_performance_files(tmp_path)
    rels = {str(p.relative_to(tmp_path)) for p in found}
    assert rels == {
        "download/norman/.abc/performance.txt",
        "download/old_ds/.def/old_ds_performance.txt",
    }


@pytest.mark.short
def test_build_record_enriches_with_lineage_and_merged_params(tmp_path):
    """A record carries stage/module/dataset plus params merged along lineage."""
    download = tmp_path / "download" / "norman" / ".68ee97fd"
    download.mkdir(parents=True)
    (download / "parameters.json").write_text(json.dumps({"uri": "http://x"}))

    split = download / "split" / "simulation" / ".a4bbadc1"
    split.mkdir(parents=True)
    (split / "parameters.json").write_text(json.dumps({"seed": 1}))

    leaf = split / "methods" / "additive" / ".default"
    perf = _write_perf(leaf)

    record = _build_record(tmp_path, perf)

    assert record["stage"] == "methods"
    assert record["module"] == "additive"
    assert record["dataset"] == "norman"
    assert record["param_hash"] == ""  # leaf is .default
    assert record["lineage"] == "download/norman/split/simulation/methods/additive"
    assert record["path"] == str(perf.relative_to(tmp_path))
    # Params from every non-default dir along the path are merged, keyed by stage.
    params = json.loads(record["params"])
    assert params == {"download": {"uri": "http://x"}, "split": {"seed": 1}}
    # Metric survived enrichment.
    assert record["s"] == 290.07


@pytest.mark.short
def test_collect_performance_end_to_end(tmp_path):
    """The command writes a combined performances.tsv consumable by dashboard."""
    from click.testing import CliRunner

    from omnibenchmark.cli.collect import collect

    _write_perf(tmp_path / "download" / "norman" / ".abc")
    _write_perf(tmp_path / "download" / "adamson" / ".def")

    result = CliRunner().invoke(collect, ["performance", "-o", str(tmp_path)])
    assert result.exit_code == 0

    out = tmp_path / "performances.tsv"
    assert out.exists()
    with open(out, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 2
    assert {r["module"] for r in rows} == {"norman", "adamson"}
    assert all(r["stage"] == "download" for r in rows)
