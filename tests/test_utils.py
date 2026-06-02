"""
Tests for utility functions.

These tests verify general utility functions used throughout the codebase.
All tests are marked as 'short' since they test simple utility functions.
"""

import pytest
from pathlib import Path
from unittest.mock import Mock

from omnibenchmark.model._merge import merge_dict_list
from omnibenchmark.core._mc_output import format_mc_output


@pytest.mark.short
class TestMergeDictList:
    """Test the merge_dict_list utility function."""

    def test_merge_simple_dicts(self):
        """Test merging simple dictionaries."""
        dict_list = [{"a": 1, "b": 2}, {"c": 3, "d": 4}]
        result = merge_dict_list(dict_list)
        expected = {"a": 1, "b": 2, "c": 3, "d": 4}
        assert result == expected

    def test_merge_overlapping_dicts_later_wins(self):
        """Test that later dictionaries override earlier ones for same keys."""
        dict_list = [{"a": 1, "b": 2}, {"b": 20, "c": 3}]
        result = merge_dict_list(dict_list)
        expected = {"a": 1, "b": 20, "c": 3}
        assert result == expected

    def test_merge_with_none_values_skipped(self):
        """Test that None values in the list are skipped."""
        dict_list = [{"a": 1}, None, {"b": 2}, None, {"c": 3}]
        result = merge_dict_list(dict_list)
        expected = {"a": 1, "b": 2, "c": 3}
        assert result == expected

    def test_merge_empty_list_returns_empty_dict(self):
        """Test that empty list returns empty dictionary."""
        result = merge_dict_list([])
        assert result == {}

    def test_merge_list_with_only_none_returns_empty_dict(self):
        """Test that list with only None values returns empty dictionary."""
        result = merge_dict_list([None, None, None])
        assert result == {}

    def test_merge_empty_dicts(self):
        """Test merging empty dictionaries."""
        dict_list = [{}, {}, {}]
        result = merge_dict_list(dict_list)
        assert result == {}

    def test_merge_mixed_value_types(self):
        """Test merging dictionaries with different value types."""
        dict_list = [
            {"str": "hello", "int": 42},
            {"list": [1, 2, 3], "bool": True},
            {"none": None, "dict": {"nested": "value"}},
        ]
        result = merge_dict_list(dict_list)
        expected = {
            "str": "hello",
            "int": 42,
            "list": [1, 2, 3],
            "bool": True,
            "none": None,
            "dict": {"nested": "value"},
        }
        assert result == expected

    def test_merge_single_dict(self):
        """Test merging a list with single dictionary."""
        dict_list = [{"a": 1, "b": 2, "c": 3}]
        result = merge_dict_list(dict_list)
        expected = {"a": 1, "b": 2, "c": 3}
        assert result == expected


@pytest.mark.short
class TestFormatMcOutput:
    """Test the format_mc_output utility function."""

    def test_format_with_path_substitution(self):
        """Test formatting when output has path with substitutions."""
        # Create a mock output object with path
        mock_output = Mock()
        mock_output.path = "metrics/{name}_results.json"
        mock_output.id = "test_collector"

        out_dir = Path("/test/output")
        collector_id = "metric_collector"

        result = format_mc_output(mock_output, out_dir, collector_id)
        expected = "/test/output/metrics/metric_collector_results.json"
        assert result == expected

    def test_format_with_input_substitution_only(self):
        """Test formatting with only {input} substitution."""
        mock_output = Mock()
        mock_output.path = "simple_output.txt"
        mock_output.id = "test_collector"

        out_dir = Path("/base/dir")
        collector_id = "collector"

        result = format_mc_output(mock_output, out_dir, collector_id)
        expected = "/base/dir/simple_output.txt"
        assert result == expected

    def test_format_no_path_uses_id(self):
        """Test formatting when output has no path, uses out_dir/id."""
        mock_output = Mock()
        mock_output.path = None
        mock_output.id = "default_collector"

        out_dir = Path("/output/directory")
        collector_id = "some_collector"

        result = format_mc_output(mock_output, out_dir, collector_id)
        expected = "/output/directory/default_collector"
        assert result == expected

    def test_format_empty_path_uses_id(self):
        """Test formatting when output has empty path, uses out_dir/id."""
        mock_output = Mock()
        mock_output.path = ""
        mock_output.id = "empty_path_collector"

        out_dir = Path("/test/dir")
        collector_id = "collector"

        result = format_mc_output(mock_output, out_dir, collector_id)
        expected = "/test/dir/empty_path_collector"
        assert result == expected

    def test_format_with_complex_substitutions(self):
        """Test formatting with multiple occurrences of substitutions."""
        mock_output = Mock()
        mock_output.path = "stage1/{name}/stage2/{name}.final"
        mock_output.id = "test_collector"

        out_dir = Path("/complex/base")
        collector_id = "multi_stage"

        result = format_mc_output(mock_output, out_dir, collector_id)
        expected = "/complex/base/stage1/multi_stage/stage2/multi_stage.final"
        assert result == expected

    def test_format_path_object_conversion(self):
        """Test that Path objects are properly converted to strings."""
        mock_output = Mock()
        mock_output.path = "converted.txt"
        mock_output.id = "path_collector"

        # Use Path object for out_dir
        out_dir = Path("/path/object/test")
        collector_id = "converter"

        result = format_mc_output(mock_output, out_dir, collector_id)
        expected = "/path/object/test/converted.txt"
        assert result == expected
        assert isinstance(result, str)
