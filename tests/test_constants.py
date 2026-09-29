"""
Tests for constants modules.

These tests verify that constants are properly defined and accessible.
All tests are marked as 'short' since they test simple constant definitions.
"""

import pytest
from omnibenchmark.constants import LayoutDesign, COMPRESSION_GZIP, DEFAULT_COMPRESSION
from omnibenchmark.core.constants import LOCAL_TIMEOUT_VAR


@pytest.mark.short
class TestLayoutDesign:
    """Test the LayoutDesign enum."""

    def test_layout_design_values(self):
        """Test that LayoutDesign enum has expected values."""
        assert LayoutDesign.Hierarchical.value == 1
        assert LayoutDesign.Spring.value == 2

    def test_layout_design_names(self):
        """Test that LayoutDesign enum has expected names."""
        assert LayoutDesign.Hierarchical.name == "Hierarchical"
        assert LayoutDesign.Spring.name == "Spring"

    def test_layout_design_iteration(self):
        """Test that LayoutDesign enum can be iterated."""
        designs = list(LayoutDesign)
        assert len(designs) == 2
        assert LayoutDesign.Hierarchical in designs
        assert LayoutDesign.Spring in designs

    def test_layout_design_string_representation(self):
        """Test string representation of LayoutDesign enum values."""
        assert str(LayoutDesign.Hierarchical) == "LayoutDesign.Hierarchical"
        assert str(LayoutDesign.Spring) == "LayoutDesign.Spring"

    def test_layout_design_comparison(self):
        """Test that LayoutDesign enum values can be compared."""
        assert LayoutDesign.Hierarchical == LayoutDesign.Hierarchical
        assert LayoutDesign.Spring == LayoutDesign.Spring
        assert LayoutDesign.Hierarchical != LayoutDesign.Spring


@pytest.mark.short
class TestBenchmarkConstants:
    """Test benchmark constants."""

    def test_local_timeout_var_defined(self):
        """Test that LOCAL_TIMEOUT_VAR is properly defined."""
        assert LOCAL_TIMEOUT_VAR == "local_task_timeout"
        assert isinstance(LOCAL_TIMEOUT_VAR, str)
        assert len(LOCAL_TIMEOUT_VAR) > 0

    def test_constants_are_final(self):
        """Test that constants behave as immutable Final types."""
        # These should not raise errors during import and usage
        timeout_var = LOCAL_TIMEOUT_VAR

        assert timeout_var is not None

        # Verify they maintain their values
        assert timeout_var == "local_task_timeout"


@pytest.mark.short
class TestCompressionConstants:
    """Test compression constants."""

    def test_compression_gzip_value(self):
        """Test that COMPRESSION_GZIP has expected value."""
        assert COMPRESSION_GZIP == "gzip"
        assert isinstance(COMPRESSION_GZIP, str)

    def test_default_compression_is_gzip(self):
        """Test that DEFAULT_COMPRESSION is set to COMPRESSION_GZIP."""
        assert DEFAULT_COMPRESSION == COMPRESSION_GZIP
        assert DEFAULT_COMPRESSION == "gzip"
