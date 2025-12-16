"""
Tests for constants modules.

These tests verify that constants are properly defined and accessible.
All tests are marked as 'short' since they test simple constant definitions.
"""

import pytest
from omnibenchmark.constants import LayoutDesign
from omnibenchmark.benchmark.constants import LOCAL_TIMEOUT_VAR, OUTPUT_PATH_PREFIX


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

    def test_output_path_prefix_defined(self):
        """Test that OUTPUT_PATH_PREFIX is properly defined."""
        expected_path = "{input}/{stage}/{module}/{params}"
        assert OUTPUT_PATH_PREFIX == expected_path
        assert isinstance(OUTPUT_PATH_PREFIX, str)

    def test_output_path_prefix_format(self):
        """Test that OUTPUT_PATH_PREFIX contains expected placeholders."""
        assert "{input}" in OUTPUT_PATH_PREFIX
        assert "{stage}" in OUTPUT_PATH_PREFIX
        assert "{module}" in OUTPUT_PATH_PREFIX
        assert "{params}" in OUTPUT_PATH_PREFIX

    def test_output_path_prefix_structure(self):
        """Test that OUTPUT_PATH_PREFIX has expected structure."""
        # Should be path-like with forward slashes
        assert "/" in OUTPUT_PATH_PREFIX
        # Should start with {input}
        assert OUTPUT_PATH_PREFIX.startswith("{input}")
        # Should have exactly 4 path components
        assert OUTPUT_PATH_PREFIX.count("/") == 3

    def test_constants_are_final(self):
        """Test that constants behave as immutable Final types."""
        # These should not raise errors during import and usage
        timeout_var = LOCAL_TIMEOUT_VAR
        path_prefix = OUTPUT_PATH_PREFIX

        assert timeout_var is not None
        assert path_prefix is not None

        # Verify they maintain their values
        assert timeout_var == "local_task_timeout"
        assert path_prefix == "{input}/{stage}/{module}/{params}"
