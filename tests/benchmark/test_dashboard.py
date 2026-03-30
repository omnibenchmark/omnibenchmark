"""Tests for dashboard export functionality."""

import pytest

from omnibenchmark.benchmark.dashboard import (
    EXPORT_AVAILABLE,
    _extract_metric_columns,
    _is_na,
    _parse_metrics,
    create_bettr_dashboard,
)


# EXPORT_AVAILABLE is always True now (no pandas required)
assert EXPORT_AVAILABLE


@pytest.fixture
def sample_performance_rows():
    """Create sample performance data as a list of dicts."""
    return [
        {
            "module": "method_a",
            "dataset": "dataset1",
            "s": 1.5,
            "max_rss": 100.0,
            "io_in": 10.0,
            "path": "/path/to/method_a/dataset1",
            "params": "default",
            "lineage": "a/b/c",
        },
        {
            "module": "method_a",
            "dataset": "dataset2",
            "s": 2.0,
            "max_rss": 120.0,
            "io_in": 15.0,
            "path": "/path/to/method_a/dataset2",
            "params": "default",
            "lineage": "a/b/c",
        },
        {
            "module": "method_b",
            "dataset": "dataset1",
            "s": 1.2,
            "max_rss": 90.0,
            "io_in": 8.0,
            "path": "/path/to/method_b/dataset1",
            "params": "default",
            "lineage": "d/e/f",
        },
        {
            "module": "method_b",
            "dataset": "dataset2",
            "s": 1.8,
            "max_rss": 110.0,
            "io_in": 12.0,
            "path": "/path/to/method_b/dataset2",
            "params": "default",
            "lineage": "d/e/f",
        },
    ]


@pytest.fixture
def minimal_performance_rows():
    """Create minimal performance data with a single metric."""
    return [
        {
            "module": "method_x",
            "dataset": "ds1",
            "s": 1.0,
            "path": "/path1",
            "params": "default",
            "lineage": "a",
        },
        {
            "module": "method_x",
            "dataset": "ds2",
            "s": 2.0,
            "path": "/path2",
            "params": "default",
            "lineage": "b",
        },
    ]


@pytest.mark.short
def test_is_na():
    """Test _is_na with various values."""
    assert _is_na(None)
    assert _is_na(float("nan"))
    assert not _is_na(0.0)
    assert not _is_na(1.5)
    assert not _is_na("value")
    assert not _is_na(0)


@pytest.mark.short
def test_extract_metric_columns(sample_performance_rows):
    """Test extraction of known performance metric columns."""
    cols = _extract_metric_columns(sample_performance_rows, id_col="module")
    assert "s" in cols
    assert "max_rss" in cols
    assert "io_in" in cols
    # Non-metric and excluded cols should not appear
    assert "module" not in cols
    assert "dataset" not in cols
    assert "path" not in cols
    assert "params" not in cols
    assert "lineage" not in cols


@pytest.mark.short
def test_extract_metric_columns_empty():
    """Test that empty input returns empty list."""
    assert _extract_metric_columns([], id_col="module") == []


@pytest.mark.short
def test_parse_metrics(sample_performance_rows):
    """Test metric parsing groups rows correctly."""
    metric_cols = ["s", "max_rss"]
    records = _parse_metrics(sample_performance_rows, metric_cols, id_col="module")

    # Two modules, each with measurement_num=1
    assert len(records) == 2
    module_ids = {r["module_id"] for r in records}
    assert "method_a_1" in module_ids
    assert "method_b_1" in module_ids

    # Each record should have metrics for both datasets
    for record in records:
        datasets = {m["dataset"] for m in record["metrics"]}
        assert "dataset1" in datasets
        assert "dataset2" in datasets


@pytest.mark.short
def test_create_bettr_dashboard_basic(sample_performance_rows):
    """Test basic bettr dashboard creation with sample data."""
    result = create_bettr_dashboard(sample_performance_rows, id_col="module")

    assert isinstance(result, dict)
    assert "idCol" in result
    assert "data" in result
    assert "idInfo" in result
    assert "metricInfo" in result
    assert "metricColors" in result
    assert "initialTransforms" in result

    assert result["idCol"] == "module"

    assert isinstance(result["data"], list)
    assert len(result["data"]) > 0

    for row in result["data"]:
        assert "module" in row

    assert isinstance(result["metricInfo"], list)
    assert len(result["metricInfo"]) > 0

    for metric in result["metricInfo"]:
        assert "Metric" in metric
        assert "Class" in metric
        assert "Dataset" in metric

    metric_names = [m["Metric"] for m in result["metricInfo"]]
    assert any("_dataset1" in name for name in metric_names)
    assert any("_dataset2" in name for name in metric_names)


@pytest.mark.short
def test_create_bettr_dashboard_minimal(minimal_performance_rows):
    """Test dashboard with minimal (single-metric) data."""
    result = create_bettr_dashboard(minimal_performance_rows, id_col="module")

    assert result["idCol"] == "module"
    assert len(result["data"]) == 1
    assert result["data"][0]["module"] == "method_x_1"

    metric_names = [m["Metric"] for m in result["metricInfo"]]
    assert "s_ds1" in metric_names
    assert "s_ds2" in metric_names


@pytest.mark.short
def test_create_bettr_dashboard_empty_raises():
    """Test that empty input raises ValueError."""
    with pytest.raises(ValueError, match="empty"):
        create_bettr_dashboard([], id_col="module")


@pytest.mark.short
def test_create_bettr_dashboard_no_metrics_raises():
    """Test that input with no known metric columns raises ValueError."""
    rows = [
        {"module": "m", "dataset": "d", "path": "/p", "params": "x", "lineage": "y"}
    ]
    with pytest.raises(ValueError, match="No metric columns"):
        create_bettr_dashboard(rows, id_col="module")


@pytest.mark.short
def test_create_bettr_dashboard_nan_values():
    """Test that NaN/None values are handled correctly."""
    rows = [
        {
            "module": "method_a",
            "dataset": "ds1",
            "s": None,
            "path": "/p",
            "params": "default",
            "lineage": "a",
        },
        {
            "module": "method_a",
            "dataset": "ds1",
            "s": float("nan"),
            "path": "/p",
            "params": "default",
            "lineage": "a",
        },
        {
            "module": "method_b",
            "dataset": "ds1",
            "s": 1.5,
            "path": "/p",
            "params": "default",
            "lineage": "b",
        },
    ]
    result = create_bettr_dashboard(rows, id_col="module")
    assert isinstance(result, dict)
    # method_b has a valid value, so it should produce data
    data_modules = [r["module"] for r in result["data"]]
    assert any("method_b" in m for m in data_modules)


@pytest.mark.short
def test_initial_transforms_populated(sample_performance_rows):
    """Test that initialTransforms contains flip/transform for known metrics."""
    result = create_bettr_dashboard(sample_performance_rows, id_col="module")
    transforms = result["initialTransforms"]
    assert len(transforms) > 0
    for col, config in transforms.items():
        assert "flip" in config
        assert "transform" in config
