"""Tests for CLI dashboard command helpers."""

import tempfile
from pathlib import Path

import pytest

from omnibenchmark.cli.dashboard import _read_performance_tsv


@pytest.mark.short
def test_read_performance_tsv_basic():
    """Read a TSV file and coerce numeric columns."""
    tsv_content = "module\tdataset\ts\tmax_rss\tpath\n"
    tsv_content += "method_a\tdataset1\t1.5\t100.0\t/path/a\n"
    tsv_content += "method_b\tdataset1\t2.3\t80.0\t/path/b\n"

    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write(tsv_content)
        tmp_path = Path(f.name)

    try:
        rows = _read_performance_tsv(tmp_path)
        assert len(rows) == 2
        assert rows[0]["module"] == "method_a"
        assert rows[0]["s"] == 1.5
        assert rows[0]["max_rss"] == 100.0
        assert rows[0]["path"] == "/path/a"
        assert rows[1]["module"] == "method_b"
        assert rows[1]["s"] == 2.3
    finally:
        tmp_path.unlink()


@pytest.mark.short
def test_read_performance_tsv_nan_values():
    """Empty/nan string values are coerced to None."""
    tsv_content = "module\tdataset\ts\n"
    tsv_content += "method_a\tds1\t\n"
    tsv_content += "method_b\tds1\tnan\n"
    tsv_content += "method_c\tds1\t1.0\n"

    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write(tsv_content)
        tmp_path = Path(f.name)

    try:
        rows = _read_performance_tsv(tmp_path)
        assert rows[0]["s"] is None  # empty string → None
        assert rows[1]["s"] is None  # "nan" → None
        assert rows[2]["s"] == 1.0
    finally:
        tmp_path.unlink()


@pytest.mark.short
def test_read_performance_tsv_string_columns():
    """Non-numeric values stay as strings."""
    tsv_content = "module\tdataset\ts\tparams\n"
    tsv_content += "method_a\tds1\t1.5\tdefault\n"

    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write(tsv_content)
        tmp_path = Path(f.name)

    try:
        rows = _read_performance_tsv(tmp_path)
        assert rows[0]["params"] == "default"  # stays as string
        assert rows[0]["s"] == 1.5  # coerced to float
    finally:
        tmp_path.unlink()
