"""
Tests for utility functions.

These tests verify general utility functions used throughout the codebase.
All tests are marked as 'short' since they test simple utility functions.
"""

import pytest
import os
from pathlib import Path
from unittest.mock import Mock, patch

from omnibenchmark.utils import (
    try_avail_envmodule,
    as_list,
    merge_dict_list,
    format_mc_output,
)


@pytest.mark.short
class TestTryAvailEnvmodule:
    """Test the try_avail_envmodule function."""

    @patch("subprocess.run")
    def test_module_available_returns_true(self, mock_run):
        """Test that available module returns True."""
        # Mock successful module check (no "No module" message in stderr)
        mock_run.return_value = Mock(
            stderr="Available Modules:\n  gcc/11.2.0\n  gcc/12.1.0"
        )

        result = try_avail_envmodule("gcc")
        assert result is True

        # Verify the command structure
        mock_run.assert_called_once()
        call_args = mock_run.call_args
        assert "module avail gcc" in call_args[0][0]
        assert call_args[1]["shell"] is True
        assert call_args[1]["text"] is True

    @patch("subprocess.run")
    def test_module_unavailable_returns_false(self, mock_run):
        """Test that unavailable module returns False."""
        # Mock failed module check
        mock_run.return_value = Mock(stderr="No module(s) or extension(s) found!")

        result = try_avail_envmodule("nonexistent-module")
        assert result is False

    @patch("subprocess.run")
    def test_partial_match_in_error_returns_false(self, mock_run):
        """Test that partial error message match still returns False."""
        mock_run.return_value = Mock(
            stderr="Some other error. No module(s) or extension(s) found! More text."
        )

        result = try_avail_envmodule("test-module")
        assert result is False

    @patch("subprocess.run")
    def test_empty_stderr_returns_true(self, mock_run):
        """Test that empty stderr (successful check) returns True."""
        mock_run.return_value = Mock(stderr="")

        result = try_avail_envmodule("some-module")
        assert result is True

    @patch("subprocess.run")
    def test_environment_variables_passed(self, mock_run):
        """Test that environment variables are properly passed to subprocess."""
        mock_run.return_value = Mock(stderr="")

        # Set a test environment variable
        test_env_var = "TEST_OMNIBENCH_VAR"
        test_env_value = "test_value"
        os.environ[test_env_var] = test_env_value

        try:
            try_avail_envmodule("test-module")

            # Check that environment was passed
            call_args = mock_run.call_args
            passed_env = call_args[1]["env"]
            assert test_env_var in passed_env
            assert passed_env[test_env_var] == test_env_value
        finally:
            # Clean up
            if test_env_var in os.environ:
                del os.environ[test_env_var]

    @patch("subprocess.run")
    def test_command_structure(self, mock_run):
        """Test that the command has the expected structure."""
        mock_run.return_value = Mock(stderr="")

        try_avail_envmodule("test-module")

        call_args = mock_run.call_args
        command = call_args[0][0]

        # Verify command components
        assert '. "$LMOD_PKG"/init/profile' in command
        assert "module purge" in command
        assert "module avail test-module" in command


@pytest.mark.short
class TestAsList:
    """Test the as_list utility function."""

    def test_list_input_returns_same_list(self):
        """Test that list input returns the same list."""
        input_list = [1, 2, 3, "test"]
        result = as_list(input_list)
        assert result == input_list
        assert result is input_list  # Should be the same object

    def test_string_input_returns_list_with_string(self):
        """Test that string input returns list containing the string."""
        input_str = "test_string"
        result = as_list(input_str)
        assert result == ["test_string"]
        assert isinstance(result, list)
        assert len(result) == 1

    def test_integer_input_returns_list_with_integer(self):
        """Test that integer input returns list containing the integer."""
        input_int = 42
        result = as_list(input_int)
        assert result == [42]
        assert isinstance(result, list)
        assert len(result) == 1

    def test_none_input_returns_list_with_none(self):
        """Test that None input returns list containing None."""
        result = as_list(None)
        assert result == [None]
        assert isinstance(result, list)
        assert len(result) == 1

    def test_dict_input_returns_list_with_dict(self):
        """Test that dict input returns list containing the dict."""
        input_dict = {"key": "value", "number": 123}
        result = as_list(input_dict)
        assert result == [input_dict]
        assert isinstance(result, list)
        assert len(result) == 1
        assert result[0] is input_dict

    def test_empty_list_returns_empty_list(self):
        """Test that empty list returns empty list."""
        input_list = []
        result = as_list(input_list)
        assert result == []
        assert result is input_list

    def test_tuple_input_returns_list_with_tuple(self):
        """Test that tuple input returns list containing the tuple."""
        input_tuple = (1, 2, 3)
        result = as_list(input_tuple)
        assert result == [(1, 2, 3)]
        assert isinstance(result, list)
        assert len(result) == 1


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
        mock_output.path = "{input}/metrics/{name}_results.json"
        mock_output.id = "test_collector"

        out_dir = Path("/test/output")
        collector_id = "metric_collector"

        result = format_mc_output(mock_output, out_dir, collector_id)
        expected = "/test/output/metrics/metric_collector_results.json"
        assert result == expected

    def test_format_with_input_substitution_only(self):
        """Test formatting with only {input} substitution."""
        mock_output = Mock()
        mock_output.path = "{input}/simple_output.txt"
        mock_output.id = "test_collector"

        out_dir = Path("/base/dir")
        collector_id = "collector"

        result = format_mc_output(mock_output, out_dir, collector_id)
        expected = "/base/dir/simple_output.txt"
        assert result == expected

    def test_format_with_name_substitution_only(self):
        """Test formatting with only {name} substitution."""
        mock_output = Mock()
        mock_output.path = "/fixed/path/{name}.log"
        mock_output.id = "test_collector"

        out_dir = Path("/test/output")
        collector_id = "log_collector"

        result = format_mc_output(mock_output, out_dir, collector_id)
        expected = "/fixed/path/log_collector.log"
        assert result == expected

    def test_format_without_substitutions(self):
        """Test formatting path without any substitution placeholders."""
        mock_output = Mock()
        mock_output.path = "/absolute/fixed/path/file.txt"
        mock_output.id = "test_collector"

        out_dir = Path("/test/output")
        collector_id = "collector"

        result = format_mc_output(mock_output, out_dir, collector_id)
        expected = "/absolute/fixed/path/file.txt"
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
        mock_output.path = "{input}/stage1/{name}/stage2/{input}/{name}.final"
        mock_output.id = "test_collector"

        out_dir = Path("/complex/base")
        collector_id = "multi_stage"

        result = format_mc_output(mock_output, out_dir, collector_id)
        expected = (
            "/complex/base/stage1/multi_stage/stage2//complex/base/multi_stage.final"
        )
        assert result == expected

    def test_format_path_object_conversion(self):
        """Test that Path objects are properly converted to strings."""
        mock_output = Mock()
        mock_output.path = "{input}/converted.txt"
        mock_output.id = "path_collector"

        # Use Path object for out_dir
        out_dir = Path("/path/object/test")
        collector_id = "converter"

        result = format_mc_output(mock_output, out_dir, collector_id)
        expected = "/path/object/test/converted.txt"
        assert result == expected
        assert isinstance(result, str)
