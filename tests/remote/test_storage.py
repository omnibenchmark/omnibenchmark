"""
Unit tests for omnibenchmark.remote.storage module.

Tests focus on behavior-driven scenarios with minimal mocking.
"""

import pytest
from unittest.mock import MagicMock

from omnibenchmark.model import Benchmark


@pytest.fixture
def mock_benchmark_with_s3():
    """Create a mock benchmark with S3 storage configuration."""
    benchmark = MagicMock(spec=Benchmark)
    benchmark.get_storage_api.return_value = "S3"
    benchmark.get_storage_endpoint.return_value = "https://s3.amazonaws.com"
    benchmark.get_storage_bucket_name.return_value = "test-bucket"
    return benchmark


@pytest.fixture
def mock_benchmark_without_storage():
    """Create a mock benchmark without storage configuration."""
    benchmark = MagicMock(spec=Benchmark)
    benchmark.get_storage_api.return_value = None
    benchmark.get_storage_endpoint.return_value = None
    benchmark.get_storage_bucket_name.return_value = None
    return benchmark


@pytest.mark.short
def test_remote_storage_args_with_credentials_in_env(
    mock_benchmark_with_s3, monkeypatch
):
    """Test remote_storage_args when credentials are available in environment."""
    from omnibenchmark.remote.storage import remote_storage_args

    # Set up environment variables
    monkeypatch.setenv("OB_STORAGE_S3_ACCESS_KEY", "test_access_key")
    monkeypatch.setenv("OB_STORAGE_S3_SECRET_KEY", "test_secret_key")

    # Call with required=False (default)
    result = remote_storage_args(mock_benchmark_with_s3, required=False)

    assert result["endpoint"] == "https://s3.amazonaws.com"
    assert result["secure"] is True
    assert result["access_key"] == "test_access_key"
    assert result["secret_key"] == "test_secret_key"


@pytest.mark.short
def test_remote_storage_args_without_credentials_not_required(
    mock_benchmark_with_s3, monkeypatch
):
    """Test remote_storage_args when credentials are missing but not required."""
    from omnibenchmark.remote.storage import remote_storage_args

    # Remove credentials from environment
    monkeypatch.delenv("OB_STORAGE_S3_ACCESS_KEY", raising=False)
    monkeypatch.delenv("OB_STORAGE_S3_SECRET_KEY", raising=False)
    monkeypatch.delenv("OB_STORAGE_S3_CONFIG", raising=False)

    # Call with required=False (should not fail)
    result = remote_storage_args(mock_benchmark_with_s3, required=False)

    # Should return base config without credentials
    assert result["endpoint"] == "https://s3.amazonaws.com"
    assert result["secure"] is True
    assert "access_key" not in result
    assert "secret_key" not in result


@pytest.mark.short
def test_remote_storage_args_without_credentials_required(
    mock_benchmark_with_s3, monkeypatch
):
    """Test remote_storage_args when credentials are missing and required."""
    from omnibenchmark.remote.storage import remote_storage_args

    # Remove credentials from environment
    monkeypatch.delenv("OB_STORAGE_S3_ACCESS_KEY", raising=False)
    monkeypatch.delenv("OB_STORAGE_S3_SECRET_KEY", raising=False)
    monkeypatch.delenv("OB_STORAGE_S3_CONFIG", raising=False)

    # Call with required=True should exit
    with pytest.raises(SystemExit) as excinfo:
        remote_storage_args(mock_benchmark_with_s3, required=True)

    assert excinfo.value.code == 1


@pytest.mark.short
def test_remote_storage_args_with_http_endpoint(mock_benchmark_with_s3, monkeypatch):
    """Test remote_storage_args with HTTP endpoint (non-secure)."""
    from omnibenchmark.remote.storage import remote_storage_args

    # Set up HTTP endpoint
    mock_benchmark_with_s3.get_storage_endpoint.return_value = "http://localhost:9000"

    # Set up credentials
    monkeypatch.setenv("OB_STORAGE_S3_ACCESS_KEY", "minioadmin")
    monkeypatch.setenv("OB_STORAGE_S3_SECRET_KEY", "minioadmin123")

    result = remote_storage_args(mock_benchmark_with_s3, required=False)

    assert result["endpoint"] == "http://localhost:9000"
    assert result["secure"] is False  # HTTP should be non-secure
    assert result["access_key"] == "minioadmin"


@pytest.mark.short
def test_remote_storage_args_with_credentials_from_file(
    mock_benchmark_with_s3, tmp_path, monkeypatch
):
    """Test remote_storage_args loading credentials from config file."""
    from omnibenchmark.remote.storage import remote_storage_args
    import json

    # Create credentials file
    creds_file = tmp_path / "s3_creds.json"
    creds = {
        "access_key": "file_access_key",
        "secret_key": "file_secret_key",
    }
    with open(creds_file, "w") as f:
        json.dump(creds, f)

    # Remove env vars and set config file
    monkeypatch.delenv("OB_STORAGE_S3_ACCESS_KEY", raising=False)
    monkeypatch.delenv("OB_STORAGE_S3_SECRET_KEY", raising=False)
    monkeypatch.setenv("OB_STORAGE_S3_CONFIG", str(creds_file))

    result = remote_storage_args(mock_benchmark_with_s3, required=False)

    assert result["access_key"] == "file_access_key"
    assert result["secret_key"] == "file_secret_key"


@pytest.mark.short
def test_remote_storage_args_with_incomplete_config_file(
    mock_benchmark_with_s3, tmp_path, monkeypatch
):
    """Test remote_storage_args with incomplete credentials in config file."""
    from omnibenchmark.remote.storage import remote_storage_args
    import json

    # Create incomplete credentials file (missing secret_key)
    creds_file = tmp_path / "s3_creds.json"
    creds = {"access_key": "file_access_key"}  # Missing secret_key
    with open(creds_file, "w") as f:
        json.dump(creds, f)

    # Set config file
    monkeypatch.delenv("OB_STORAGE_S3_ACCESS_KEY", raising=False)
    monkeypatch.delenv("OB_STORAGE_S3_SECRET_KEY", raising=False)
    monkeypatch.setenv("OB_STORAGE_S3_CONFIG", str(creds_file))

    # With required=False, should return empty credentials
    result = remote_storage_args(mock_benchmark_with_s3, required=False)
    assert "access_key" not in result or "secret_key" not in result

    # With required=True, should exit
    with pytest.raises(SystemExit):
        remote_storage_args(mock_benchmark_with_s3, required=True)


@pytest.mark.short
def test_remote_storage_args_no_storage_configured(mock_benchmark_without_storage):
    """Test remote_storage_args when benchmark has no storage configured."""
    from omnibenchmark.remote.storage import remote_storage_args

    result = remote_storage_args(mock_benchmark_without_storage, required=False)

    assert result == {}


@pytest.mark.short
def test_remote_storage_args_minio_api(mock_benchmark_with_s3, monkeypatch):
    """Test remote_storage_args with MinIO API (should work same as S3)."""
    from omnibenchmark.remote.storage import remote_storage_args

    # Set API to MinIO instead of S3
    mock_benchmark_with_s3.get_storage_api.return_value = "MinIO"

    monkeypatch.setenv("OB_STORAGE_S3_ACCESS_KEY", "test_key")
    monkeypatch.setenv("OB_STORAGE_S3_SECRET_KEY", "test_secret")

    result = remote_storage_args(mock_benchmark_with_s3, required=False)

    # Should work the same as S3
    assert result["access_key"] == "test_key"
    assert result["secret_key"] == "test_secret"


@pytest.mark.short
def test_remote_storage_args_missing_endpoint(mock_benchmark_with_s3, monkeypatch):
    """Test remote_storage_args when endpoint is None."""
    from omnibenchmark.remote.storage import remote_storage_args

    # Set endpoint to None
    mock_benchmark_with_s3.get_storage_endpoint.return_value = None

    monkeypatch.setenv("OB_STORAGE_S3_ACCESS_KEY", "test_key")
    monkeypatch.setenv("OB_STORAGE_S3_SECRET_KEY", "test_secret")

    result = remote_storage_args(mock_benchmark_with_s3, required=False)

    # Should return empty dict when endpoint is missing
    assert result == {}


@pytest.mark.short
def test_remote_storage_args_required_parameter_enforces_validation(
    mock_benchmark_with_s3, monkeypatch
):
    """Test that required=True parameter actually enforces credential validation."""
    from omnibenchmark.remote.storage import remote_storage_args

    # Remove all credentials
    monkeypatch.delenv("OB_STORAGE_S3_ACCESS_KEY", raising=False)
    monkeypatch.delenv("OB_STORAGE_S3_SECRET_KEY", raising=False)
    monkeypatch.delenv("OB_STORAGE_S3_CONFIG", raising=False)

    # With required=False, should return config without credentials
    result = remote_storage_args(mock_benchmark_with_s3, required=False)
    assert result["endpoint"] == "https://s3.amazonaws.com"
    assert "access_key" not in result or not result.get("access_key")

    # With required=True, should exit with error
    with pytest.raises(SystemExit) as excinfo:
        remote_storage_args(mock_benchmark_with_s3, required=True)
    assert excinfo.value.code == 1


@pytest.mark.short
def test_remote_storage_args_respects_required_flag_with_credentials(
    mock_benchmark_with_s3, monkeypatch
):
    """Test that required flag works correctly when credentials ARE present."""
    from omnibenchmark.remote.storage import remote_storage_args

    # Set up credentials
    monkeypatch.setenv("OB_STORAGE_S3_ACCESS_KEY", "test_key")
    monkeypatch.setenv("OB_STORAGE_S3_SECRET_KEY", "test_secret")

    # Both required=True and required=False should work with credentials present
    result_not_required = remote_storage_args(mock_benchmark_with_s3, required=False)
    result_required = remote_storage_args(mock_benchmark_with_s3, required=True)

    # Both should return the same result with credentials
    assert result_not_required["access_key"] == "test_key"
    assert result_required["access_key"] == "test_key"
    assert result_not_required == result_required
