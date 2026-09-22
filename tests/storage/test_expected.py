"""The extra files archiving keeps beyond the benchmark's declared outputs."""

import pytest

from omnibenchmark.storage import StorageOptions, get_expected_benchmark_output_files


class _FakeBenchmark:
    def get_output_paths(self):
        return set()


@pytest.mark.short
def test_fan_in_lineage_sidecar_is_archived(tmp_path, monkeypatch):
    """lineage.json is the only record of a fan-in node's ancestry (010 §5.2).

    It is not a declared Snakemake output, so if the glob list drops it an
    archived gather result can never say what it was computed from.
    """
    node = tmp_path / "out" / "metrics" / "D1" / "MC" / ".default"
    node.mkdir(parents=True)
    (node / "lineage.json").write_text("{}")
    (node / "parameters.json").write_text("{}")
    monkeypatch.chdir(tmp_path)

    found = get_expected_benchmark_output_files(_FakeBenchmark(), StorageOptions("out"))

    assert "out/metrics/D1/MC/.default/lineage.json" in found
    assert "out/metrics/D1/MC/.default/parameters.json" in found
