"""Tests for remote policy create command."""

import json
from click.testing import CliRunner
from omnibenchmark.cli.remote import create_policy


class TestCreatePolicyCommand:
    """Tests for policy creation command."""

    def test_create_policy_with_bucket_name(self):
        """Test creating policy with direct bucket name (no YAML file needed)."""
        runner = CliRunner()
        result = runner.invoke(create_policy, ["--bucket", "test-bucket"])

        assert result.exit_code == 0

        # Parse the output as JSON
        policy = json.loads(result.output)

        # Verify policy structure
        assert policy["Version"] == "2012-10-17"
        assert len(policy["Statement"]) == 2

        # Verify bucket name in resources
        assert "arn:aws:s3:::test-bucket/*" in policy["Statement"][0]["Resource"]
        assert "arn:aws:s3:::test-bucket" in policy["Statement"][0]["Resource"]

        # Verify allow statement
        assert policy["Statement"][0]["Effect"] == "Allow"
        assert "s3:*" in policy["Statement"][0]["Action"]

        # Verify deny statement
        assert policy["Statement"][1]["Effect"] == "Deny"
        assert "s3:BypassGovernanceRetention" in policy["Statement"][1]["Action"]

    def test_create_policy_from_yaml(self, tmp_path):
        """Test creating policy from benchmark YAML file."""
        # Create a test YAML file
        yaml_content = """id: test_benchmark
description: Test
version: "1.0"
benchmarker: "Test"
storage:
  api: S3
  endpoint: https://example.com
  bucket_name: my-yaml-bucket
software_backend: host
software_environments:
  - id: test
    description: test
stages:
  - id: data
    modules:
      - id: D1
        software_environment: test
        repository:
          url: https://example.com/repo.git
          commit: abc123
    outputs:
      - id: data.out
        path: "{dataset}.txt"
"""
        yaml_file = tmp_path / "test.yaml"
        yaml_file.write_text(yaml_content)

        runner = CliRunner()
        result = runner.invoke(create_policy, ["--benchmark", str(yaml_file)])

        assert result.exit_code == 0

        # Parse the output as JSON
        policy = json.loads(result.output)

        # Verify bucket name from YAML is used
        assert "arn:aws:s3:::my-yaml-bucket/*" in policy["Statement"][0]["Resource"]
        assert "arn:aws:s3:::my-yaml-bucket" in policy["Statement"][0]["Resource"]

    def test_create_policy_yaml_without_bucket(self, tmp_path):
        """Test error when YAML has no bucket configured."""
        yaml_content = """id: test_benchmark
description: Test
version: "1.0"
benchmarker: "Test"
storage:
  api: S3
  endpoint: https://example.com
software_backend: host
software_environments:
  - id: test
    description: test
stages:
  - id: data
    modules:
      - id: D1
        software_environment: test
        repository:
          url: https://example.com/repo.git
          commit: abc123
    outputs:
      - id: data.out
        path: "{dataset}.txt"
"""
        yaml_file = tmp_path / "test.yaml"
        yaml_file.write_text(yaml_content)

        runner = CliRunner()
        result = runner.invoke(create_policy, ["--benchmark", str(yaml_file)])

        # Should exit with error code when bucket is not configured
        assert result.exit_code == 1

    def test_create_policy_no_arguments(self):
        """Test error when neither --bucket nor --benchmark provided."""
        runner = CliRunner()
        result = runner.invoke(create_policy, [])

        # Should exit with error code when no arguments provided
        assert result.exit_code == 1
