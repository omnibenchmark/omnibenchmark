"""
S3 configuration validation tests that don't require actual S3 setup.

These tests validate the S3 configuration files and basic CLI behavior
without needing MinIO or remote S3 credentials.
"""

import pytest
import yaml
from pathlib import Path


@pytest.fixture
def s3_config_path():
    """Get the path to the S3 config."""
    return Path(__file__).parent / "configs" / "06_s3_remote.yaml"


def test_s3_config_file_exists(s3_config_path):
    """Test that the S3 config file exists and is readable."""
    assert s3_config_path.exists(), f"S3 config file not found: {s3_config_path}"
    assert s3_config_path.is_file(), f"S3 config path is not a file: {s3_config_path}"


def test_s3_config_yaml_valid(s3_config_path):
    """Test that the S3 config file contains valid YAML."""
    with open(s3_config_path, "r") as f:
        config = yaml.safe_load(f)

    assert isinstance(config, dict), "Config file should contain a YAML dictionary"
    assert config, "Config file should not be empty"


def test_s3_config_has_required_fields(s3_config_path):
    """Test that the S3 config has all required fields."""
    with open(s3_config_path, "r") as f:
        config = yaml.safe_load(f)

    # Required top-level fields
    required_fields = ["id", "description", "version", "storage", "stages"]
    for field in required_fields:
        assert field in config, f"Missing required field: {field}"

    # Storage configuration validation
    storage = config["storage"]
    assert storage["api"] == "S3", "Storage API should be S3"
    assert "bucket_name" in storage, "Storage should have bucket_name"
    assert "endpoint" in storage, "Storage should have endpoint"

    # Stages validation
    stages = config["stages"]
    assert len(stages) >= 1, "Should have at least one stage"

    # Check for data stage
    data_stage = next((s for s in stages if s["id"] == "data"), None)
    assert data_stage is not None, "Should have a data stage"
    assert "modules" in data_stage, "Data stage should have modules"
    assert len(data_stage["modules"]) >= 1, "Data stage should have at least one module"


def test_s3_config_minimal_structure(s3_config_path):
    """Test that the S3 config has the expected minimal structure for fast testing."""
    with open(s3_config_path, "r") as f:
        config = yaml.safe_load(f)

    # Should have reasonable structure for testing
    stages = config["stages"]
    assert len(stages) >= 1, "Should have at least 1 stage"
    assert len(stages) <= 5, "Should have at most 5 stages for reasonable testing"

    # Data stage should exist
    data_stage = next((s for s in stages if s["id"] == "data"), None)
    assert data_stage is not None, "Should have a data stage"
    data_modules = data_stage["modules"]
    assert len(data_modules) >= 1, "Data stage should have at least 1 module"

    # Check parameter values are reasonable for testing
    for module in data_modules:
        if "parameters" in module:
            for param in module["parameters"]:
                if "evaluate" in param and isinstance(param["evaluate"], str):
                    try:
                        value = int(param["evaluate"])
                        assert (
                            value <= 1000
                        ), f"Parameter values should be reasonable for testing, got {value}"
                    except ValueError:
                        pass  # Skip non-numeric parameters


def test_s3_config_software_environment(s3_config_path):
    """Test that the S3 config uses host environment for simplicity."""
    with open(s3_config_path, "r") as f:
        config = yaml.safe_load(f)

    # Should use host environment
    assert "software_environments" in config, "Should have software_environments"
    assert "host" in config["software_environments"], "Should have host environment"

    # All modules should use host environment
    for stage in config["stages"]:
        for module in stage["modules"]:
            assert (
                module["software_environment"] == "host"
            ), "All modules should use host environment for simplicity"


def test_s3_config_bucket_name_format(s3_config_path):
    """Test that the S3 config has appropriate bucket name format."""
    with open(s3_config_path, "r") as f:
        config = yaml.safe_load(f)

    bucket_name = config["storage"]["bucket_name"]

    # Should be a reasonable bucket name
    assert isinstance(bucket_name, str), "Bucket name should be a string"
    assert len(bucket_name) >= 3, "Bucket name should be at least 3 characters"
    assert len(bucket_name) <= 63, "Bucket name should be at most 63 characters"
    assert (
        bucket_name.islower() or "-" in bucket_name
    ), "Bucket name should follow S3 naming conventions"

    # Bucket name can be generic (actual bucket name is generated with timestamp at runtime)
    # Just verify it follows basic S3 naming rules
    assert (
        bucket_name.replace("-", "").replace("_", "").isalnum()
        or bucket_name.replace("-", "").isalnum()
    ), "Bucket name should only contain alphanumeric characters and hyphens"


def test_s3_config_endpoint_format(s3_config_path):
    """Test that the S3 config has valid endpoint format."""
    with open(s3_config_path, "r") as f:
        config = yaml.safe_load(f)

    endpoint = config["storage"]["endpoint"]

    # Should be a valid URL-like string
    assert isinstance(endpoint, str), "Endpoint should be a string"
    assert endpoint.startswith("http://") or endpoint.startswith(
        "https://"
    ), "Endpoint should be a valid URL"

    # For local testing, should use localhost
    if "localhost" in endpoint:
        assert ":9000" in endpoint, "Local endpoint should use port 9000 for MinIO"


@pytest.mark.skip(reason="Expected results file not used for S3 tests")
def test_expected_results_file_exists():
    """Test that the expected results file exists for validation."""
    expected_results_path = (
        Path(__file__).parent / "configs" / "06_s3_remote.expected.json"
    )
    assert (
        expected_results_path.exists()
    ), f"Expected results file not found: {expected_results_path}"

    import json

    with open(expected_results_path, "r") as f:
        expected = json.load(f)

    assert "test_name" in expected, "Expected results should have test_name"
    assert (
        expected["test_name"] == "06_s3_minimal"
    ), "Expected results should match test name"
    assert (
        "expected_files" in expected
    ), "Expected results should have expected_files list"
    assert (
        len(expected["expected_files"]) >= 1
    ), "Should expect at least one output file"


def test_docker_compose_file_exists():
    """Test that the docker-compose setup files exist."""
    s3_dir = Path(__file__).parent / "s3"

    docker_compose_path = s3_dir / "docker-compose.yml"
    assert (
        docker_compose_path.exists()
    ), f"Docker compose file not found: {docker_compose_path}"

    setup_script_path = s3_dir / "setup-minio.sh"
    assert (
        setup_script_path.exists()
    ), f"MinIO setup script not found: {setup_script_path}"

    run_script_path = s3_dir / "run-local-s3-test.sh"
    assert (
        run_script_path.exists()
    ), f"Local test runner script not found: {run_script_path}"


def test_readme_documentation_exists():
    """Test that documentation exists for the S3 setup."""
    s3_dir = Path(__file__).parent / "s3"
    readme_path = s3_dir / "README.md"

    assert readme_path.exists(), f"README file not found: {readme_path}"

    # Check that README contains key sections
    with open(readme_path, "r") as f:
        content = f.read()

    required_sections = [
        "Quick Start",
        "Local Testing with MinIO",
        "Remote S3 Testing",
        "Environment Variables",
        "Troubleshooting",
    ]

    for section in required_sections:
        assert section in content, f"README should contain section: {section}"
