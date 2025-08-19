"""Tests for strict semantic version validation in Benchmark model."""

import pytest
from pydantic import ValidationError

from omnibenchmark.model import Benchmark


class TestVersionValidation:
    """Test suite for version field validation."""

    @pytest.mark.short
    def test_valid_three_part_versions(self):
        """Test valid x.y.z format versions."""
        valid_versions = [
            "0.0.0",
            "0.0.1",
            "0.1.0",
            "1.0.0",
            "1.2.3",
            "10.20.30",
            "999.999.999",
            "2.0.0",
            "3.14.159",
        ]

        for version in valid_versions:
            benchmark_data = self._get_minimal_benchmark_data()
            benchmark_data["version"] = version
            benchmark = Benchmark(**benchmark_data)
            assert benchmark.version == version

    @pytest.mark.short
    def test_valid_two_part_versions(self):
        """Test valid x.y format versions."""
        valid_versions = [
            "0.0",
            "0.1",
            "1.0",
            "1.2",
            "10.20",
            "999.999",
            "2.0",
            "3.14",
        ]

        for version in valid_versions:
            benchmark_data = self._get_minimal_benchmark_data()
            benchmark_data["version"] = version
            benchmark = Benchmark(**benchmark_data)
            assert benchmark.version == version

    @pytest.mark.short
    def test_numeric_version_rejected(self):
        """Test that numeric values are rejected (no auto-conversion)."""
        benchmark_data = self._get_minimal_benchmark_data()

        # Float version should be rejected
        benchmark_data["version"] = 1.0
        with pytest.raises(ValidationError) as exc_info:
            Benchmark(**benchmark_data)
        assert "Input should be a valid string" in str(exc_info.value)

        # Integer version should be rejected
        benchmark_data["version"] = 2
        with pytest.raises(ValidationError) as exc_info:
            Benchmark(**benchmark_data)
        assert "Input should be a valid string" in str(exc_info.value)

    @pytest.mark.short
    def test_invalid_version_formats(self):
        """Test that invalid version formats are rejected."""
        invalid_versions = [
            "1",  # Single number
            "1.",  # Trailing dot
            ".1",  # Leading dot
            "1.2.",  # Trailing dot after two parts
            "1.2.3.4",  # Four parts
            "1.2.3-alpha",  # Pre-release version
            "1.2.3+build",  # Build metadata
            "1.2.3-alpha+build",  # Pre-release and build
            "v1.2.3",  # Version prefix
            "1.2.a",  # Non-numeric patch
            "a.b.c",  # All non-numeric
            "1,2,3",  # Wrong separator
            "1_2_3",  # Underscore separator
            "-1.2.3",  # Negative major
            "1.-2.3",  # Negative minor
            "1.2.-3",  # Negative patch
            "01.2.3",  # Leading zero in major
            "1.02.3",  # Leading zero in minor
            "1.2.03",  # Leading zero in patch
            "",  # Empty string
            "latest",  # Text version
        ]

        for version in invalid_versions:
            benchmark_data = self._get_minimal_benchmark_data()
            benchmark_data["version"] = version
            with pytest.raises(ValidationError) as exc_info:
                Benchmark(**benchmark_data)
            # Check for either semantic versioning error or string type error
            error_msg = str(exc_info.value)
            assert (
                "does not follow strict semantic versioning" in error_msg
                or "Input should be a valid string" in error_msg
                or "must be a non-empty string" in error_msg
            )

    @pytest.mark.short
    def test_single_integer_rejected(self):
        """Test that single integers are rejected (no conversion)."""
        benchmark_data = self._get_minimal_benchmark_data()
        benchmark_data["version"] = 1  # Integer should be rejected

        with pytest.raises(ValidationError) as exc_info:
            Benchmark(**benchmark_data)
        assert "Input should be a valid string" in str(exc_info.value)

    @pytest.mark.short
    def test_benchmark_yaml_spec_string_only(self):
        """Test that benchmark_yaml_spec only accepts strings (no numeric coercion)."""
        benchmark_data = self._get_minimal_benchmark_data()

        # Float yaml spec should be rejected
        benchmark_data["benchmark_yaml_spec"] = 0.01
        with pytest.raises(ValidationError) as exc_info:
            Benchmark(**benchmark_data)
        assert "Input should be a valid string" in str(exc_info.value)

        # Integer yaml spec should be rejected
        benchmark_data["benchmark_yaml_spec"] = 1
        with pytest.raises(ValidationError) as exc_info:
            Benchmark(**benchmark_data)
        assert "Input should be a valid string" in str(exc_info.value)

        # Valid string yaml spec should work
        benchmark_data["benchmark_yaml_spec"] = "0.01"
        benchmark = Benchmark(**benchmark_data)
        assert benchmark.benchmark_yaml_spec == "0.01"

        # None is allowed
        benchmark_data["benchmark_yaml_spec"] = None
        benchmark = Benchmark(**benchmark_data)
        assert benchmark.benchmark_yaml_spec is None

    @pytest.mark.short
    def test_version_with_leading_zeros_rejected(self):
        """Test that versions with leading zeros are rejected."""
        invalid_versions_with_leading_zeros = [
            "01.0.0",
            "1.01.0",
            "1.0.01",
            "01.01",
            "1.01",
        ]

        for version in invalid_versions_with_leading_zeros:
            benchmark_data = self._get_minimal_benchmark_data()
            benchmark_data["version"] = version
            with pytest.raises(ValidationError) as exc_info:
                Benchmark(**benchmark_data)
            error_msg = str(exc_info.value)
            assert (
                "does not follow strict semantic versioning" in error_msg
                or "must be a non-empty string" in error_msg
            )

    @pytest.mark.short
    def test_zero_versions_allowed(self):
        """Test that versions with zero as major/minor/patch are allowed."""
        valid_zero_versions = [
            "0.0.0",
            "0.1.0",
            "1.0.0",
            "1.1.0",
            "0.0",
            "1.0",
        ]

        for version in valid_zero_versions:
            benchmark_data = self._get_minimal_benchmark_data()
            benchmark_data["version"] = version
            benchmark = Benchmark(**benchmark_data)
            assert benchmark.version == version

    def _get_minimal_benchmark_data(self):
        """Helper to get minimal valid benchmark data for testing."""
        return {
            "id": "test_benchmark",
            "name": "Test Benchmark",
            "benchmarker": "Test Author",
            "version": "1.0.0",  # Default valid version
            "software_backend": "host",
            "software_environments": [],
            "stages": [],
        }
