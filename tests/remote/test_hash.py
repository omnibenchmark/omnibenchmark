"""
Tests for hash utility functions.

These tests verify hash computation functionality for files.
All tests are marked as 'short' since they test simple utility functions.
"""

import pytest
import tempfile
import hashlib
from pathlib import Path
from omnibenchmark.remote.hash import checksum


@pytest.mark.short
class TestChecksum:
    """Test the checksum function."""

    def test_checksum_empty_file(self):
        """Test checksum of empty file."""
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as f:
            temp_path = f.name

        try:
            result = checksum(temp_path)
            # MD5 of empty file is always d41d8cd98f00b204e9800998ecf8427e
            expected = "d41d8cd98f00b204e9800998ecf8427e"
            assert result == expected
        finally:
            Path(temp_path).unlink()

    def test_checksum_simple_content(self):
        """Test checksum of file with simple content."""
        content = "hello world"

        with tempfile.NamedTemporaryFile(mode="w", delete=False) as f:
            f.write(content)
            temp_path = f.name

        try:
            result = checksum(temp_path)
            # Calculate expected MD5
            expected = hashlib.md5(content.encode()).hexdigest()
            assert result == expected
        finally:
            Path(temp_path).unlink()

    def test_checksum_binary_content(self):
        """Test checksum of file with binary content."""
        binary_content = b"\x00\x01\x02\x03\xff\xfe\xfd"

        with tempfile.NamedTemporaryFile(mode="wb", delete=False) as f:
            f.write(binary_content)
            temp_path = f.name

        try:
            result = checksum(temp_path)
            expected = hashlib.md5(binary_content).hexdigest()
            assert result == expected
        finally:
            Path(temp_path).unlink()

    def test_checksum_large_file(self):
        """Test checksum of larger file to verify chunked reading."""
        # Create content larger than 4096 bytes to test chunked reading
        content = "A" * 8192  # 8KB

        with tempfile.NamedTemporaryFile(mode="w", delete=False) as f:
            f.write(content)
            temp_path = f.name

        try:
            result = checksum(temp_path)
            expected = hashlib.md5(content.encode()).hexdigest()
            assert result == expected
        finally:
            Path(temp_path).unlink()

    def test_checksum_multiline_content(self):
        """Test checksum of file with multiline content."""
        content = "line 1\nline 2\nline 3\n"

        with tempfile.NamedTemporaryFile(mode="w", delete=False) as f:
            f.write(content)
            temp_path = f.name

        try:
            result = checksum(temp_path)
            expected = hashlib.md5(content.encode()).hexdigest()
            assert result == expected
        finally:
            Path(temp_path).unlink()

    def test_checksum_unicode_content(self):
        """Test checksum of file with unicode content."""
        content = "Hello 世界! 🌍"

        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", delete=False) as f:
            f.write(content)
            temp_path = f.name

        try:
            result = checksum(temp_path)
            expected = hashlib.md5(content.encode("utf-8")).hexdigest()
            assert result == expected
        finally:
            Path(temp_path).unlink()

    def test_checksum_returns_hexdigest(self):
        """Test that checksum returns a proper hex digest."""
        content = "test content"

        with tempfile.NamedTemporaryFile(mode="w", delete=False) as f:
            f.write(content)
            temp_path = f.name

        try:
            result = checksum(temp_path)
            # Should be 32 character hex string for MD5
            assert len(result) == 32
            assert all(c in "0123456789abcdef" for c in result)
        finally:
            Path(temp_path).unlink()

    def test_checksum_consistent_results(self):
        """Test that checksum returns consistent results for same content."""
        content = "consistent content test"

        # Create same content in two different files
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as f1:
            f1.write(content)
            temp_path1 = f1.name

        with tempfile.NamedTemporaryFile(mode="w", delete=False) as f2:
            f2.write(content)
            temp_path2 = f2.name

        try:
            result1 = checksum(temp_path1)
            result2 = checksum(temp_path2)
            assert result1 == result2
        finally:
            Path(temp_path1).unlink()
            Path(temp_path2).unlink()

    def test_checksum_different_content_different_hash(self):
        """Test that different content produces different checksums."""
        content1 = "content one"
        content2 = "content two"

        with tempfile.NamedTemporaryFile(mode="w", delete=False) as f1:
            f1.write(content1)
            temp_path1 = f1.name

        with tempfile.NamedTemporaryFile(mode="w", delete=False) as f2:
            f2.write(content2)
            temp_path2 = f2.name

        try:
            result1 = checksum(temp_path1)
            result2 = checksum(temp_path2)
            assert result1 != result2
        finally:
            Path(temp_path1).unlink()
            Path(temp_path2).unlink()

    def test_checksum_file_not_found(self):
        """Test that checksum raises appropriate error for non-existent file."""
        with pytest.raises(FileNotFoundError):
            checksum("/nonexistent/file/path.txt")
