"""Unit tests for temporary directory creation in cli/run module."""

import os
import re
import getpass
from datetime import datetime

import pytest

from omnibenchmark.cli.run import create_temp_dir_with_prefix


@pytest.mark.short
class TestCreateTempDirWithPrefix:
    """Test the create_temp_dir_with_prefix function."""

    def test_creates_directory(self):
        """Test that function creates a directory."""
        temp_dir = create_temp_dir_with_prefix()
        try:
            assert os.path.exists(temp_dir), "Directory should exist"
            assert os.path.isdir(temp_dir), "Path should be a directory"
        finally:
            if os.path.exists(temp_dir):
                os.rmdir(temp_dir)

    def test_prefix_format(self):
        """Test that directory name has correct prefix format."""
        temp_dir = create_temp_dir_with_prefix()
        try:
            basename = os.path.basename(temp_dir)
            assert basename.startswith(
                "ob-run-"
            ), f"Expected prefix 'ob-run-' but got: {basename}"
        finally:
            if os.path.exists(temp_dir):
                os.rmdir(temp_dir)

    def test_contains_username(self):
        """Test that directory name contains the current username."""
        temp_dir = create_temp_dir_with_prefix()
        try:
            basename = os.path.basename(temp_dir)
            username = getpass.getuser()
            assert (
                username in basename
            ), f"Expected username '{username}' in directory name: {basename}"
        finally:
            if os.path.exists(temp_dir):
                os.rmdir(temp_dir)

    def test_contains_timestamp(self):
        """Test that directory name contains a timestamp."""
        temp_dir = create_temp_dir_with_prefix()
        try:
            basename = os.path.basename(temp_dir)
            # Check for timestamp pattern YYYYMMDD-HHMMSS
            timestamp_pattern = r"\d{8}-\d{6}"
            assert re.search(
                timestamp_pattern, basename
            ), f"Expected timestamp pattern in: {basename}"
        finally:
            if os.path.exists(temp_dir):
                os.rmdir(temp_dir)

    def test_full_prefix_format(self):
        """Test the complete prefix format: ob-run-{user}-{timestamp}-"""
        temp_dir = create_temp_dir_with_prefix()
        try:
            basename = os.path.basename(temp_dir)
            username = getpass.getuser()

            # Expected format: ob-run-{user}-{timestamp}-{random}
            # Pattern: ob-run-username-YYYYMMDD-HHMMSS-randomchars
            pattern = rf"^ob-run-{re.escape(username)}-\d{{8}}-\d{{6}}-\w+$"
            assert re.match(
                pattern, basename
            ), f"Directory name '{basename}' does not match expected pattern"
        finally:
            if os.path.exists(temp_dir):
                os.rmdir(temp_dir)

    def test_timestamp_is_current(self):
        """Test that the timestamp in the directory name is close to current time."""
        before = datetime.now()
        temp_dir = create_temp_dir_with_prefix()
        after = datetime.now()

        try:
            basename = os.path.basename(temp_dir)

            # Extract timestamp from basename
            # Format: ob-run-{user}-YYYYMMDD-HHMMSS-{random}
            timestamp_match = re.search(r"(\d{8})-(\d{6})", basename)
            assert timestamp_match, f"Could not find timestamp in: {basename}"

            date_str, time_str = timestamp_match.groups()
            timestamp_str = f"{date_str}-{time_str}"
            timestamp = datetime.strptime(timestamp_str, "%Y%m%d-%H%M%S")

            # Verify timestamp is between before and after (with 1 second tolerance)
            assert (
                before.replace(microsecond=0)
                <= timestamp
                <= after.replace(microsecond=0) + datetime.resolution
            ), f"Timestamp {timestamp} not within expected range [{before}, {after}]"
        finally:
            if os.path.exists(temp_dir):
                os.rmdir(temp_dir)

    def test_creates_unique_directories(self):
        """Test that multiple calls create unique directories."""
        temp_dirs = []
        try:
            # Create multiple temp directories
            for _ in range(3):
                temp_dir = create_temp_dir_with_prefix()
                temp_dirs.append(temp_dir)

            # Verify all directories exist
            for temp_dir in temp_dirs:
                assert os.path.exists(temp_dir), f"Directory {temp_dir} should exist"

            # Verify all directories are unique
            assert len(set(temp_dirs)) == len(
                temp_dirs
            ), "All directories should be unique"
        finally:
            # Clean up
            for temp_dir in temp_dirs:
                if os.path.exists(temp_dir):
                    os.rmdir(temp_dir)

    def test_created_in_tmp_directory(self):
        """Test that directory is created in the system's temp directory."""
        temp_dir = create_temp_dir_with_prefix()
        try:
            import tempfile

            system_tmp = tempfile.gettempdir()
            assert temp_dir.startswith(
                system_tmp
            ), f"Directory {temp_dir} should be in system temp directory {system_tmp}"
        finally:
            if os.path.exists(temp_dir):
                os.rmdir(temp_dir)
