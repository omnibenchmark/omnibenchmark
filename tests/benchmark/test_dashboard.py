"""Tests for dashboard export functionality."""

import pytest

from omnibenchmark.benchmark.dashboard import (
    create_bettr_dashboard,
    EXPORT_AVAILABLE,
)


pytestmark = pytest.mark.skipif(
    not EXPORT_AVAILABLE, reason="export not installed (export optional dependency)"
)


@pytest.fixture
def sample_performance_df():
    import pandas as pd

    """Create a sample performance DataFrame for testing."""
    data = {
        "module": ["method_a", "method_a", "method_b", "method_b"],
        "dataset": ["dataset1", "dataset2", "dataset1", "dataset2"],
        "s": [1.5, 2.0, 1.2, 1.8],
        "max_rss": [100.0, 120.0, 90.0, 110.0],
        "io_in": [10.0, 15.0, 8.0, 12.0],
        "path": [
            "/path/to/method_a/dataset1",
            "/path/to/method_a/dataset2",
            "/path/to/method_b/dataset1",
            "/path/to/method_b/dataset2",
        ],
        "params": ["default", "default", "default", "default"],
        "lineage": ["a/b/c", "a/b/c", "d/e/f", "d/e/f"],
    }
    return pd.DataFrame(data)


@pytest.fixture
def minimal_performance_df():
    import pandas as pd

    """Create a minimal performance DataFrame with single metric."""
    data = {
        "module": ["method_x", "method_x"],
        "dataset": ["ds1", "ds2"],
        "s": [1.0, 2.0],
        "path": ["/path1", "/path2"],
        "params": ["default", "default"],
        "lineage": ["a", "b"],
    }
    return pd.DataFrame(data)


@pytest.mark.short
def test_create_bettr_dashboard_basic(sample_performance_df):
    """Test basic bettr dashboard creation with sample data."""

    result = create_bettr_dashboard(sample_performance_df, id_col="module")

    # Verify basic structure
    assert isinstance(result, dict)
    assert "idCol" in result
    assert "data" in result
    assert "idInfo" in result
    assert "metricInfo" in result
    assert "metricColors" in result
    assert "initialTransforms" in result

    # Verify idCol
    assert result["idCol"] == "module"

    # Verify data structure
    assert isinstance(result["data"], list)
    assert len(result["data"]) > 0

    # Each data row should have the idCol
    for row in result["data"]:
        assert "module" in row

    # Verify metricInfo
    assert isinstance(result["metricInfo"], list)
    assert len(result["metricInfo"]) > 0

    # Each metric should have expected fields
    for metric in result["metricInfo"]:
        assert "Metric" in metric
        assert "Class" in metric
        assert "Dataset" in metric

    # Verify metric_dataset naming convention
    metric_names = [m["Metric"] for m in result["metricInfo"]]
    assert any("_dataset1" in name for name in metric_names)
    assert any("_dataset2" in name for name in metric_names)
