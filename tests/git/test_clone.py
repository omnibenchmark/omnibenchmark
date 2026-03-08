"""Short unit tests for git/clone.py pure functions and null streams."""

import pytest

from omnibenchmark.git.clone import (
    is_commit_hash,
    _NullBinaryStream,
)


@pytest.mark.short
class TestIsCommitHash:
    def test_full_sha(self):
        assert is_commit_hash("a" * 40)

    def test_short_sha_7(self):
        assert is_commit_hash("d807c924")  # 8 hex chars

    def test_short_sha_minimum(self):
        assert is_commit_hash("abcdef1")  # exactly 7

    def test_branch_name(self):
        assert not is_commit_hash("main")

    def test_tag_name(self):
        assert not is_commit_hash("v1.0.0")

    def test_too_short(self):
        assert not is_commit_hash("abc12")  # 5 chars — below minimum

    def test_uppercase_accepted(self):
        assert is_commit_hash("ABCDEF1234567")

    def test_mixed_case(self):
        assert is_commit_hash("AbCdEf1234567")

    def test_non_hex_char(self):
        assert not is_commit_hash("xyz1234567890")

    def test_empty(self):
        assert not is_commit_hash("")


@pytest.mark.short
class TestNullBinaryStream:
    def test_write_returns_length(self):
        s = _NullBinaryStream()
        assert s.write(b"hello") == 5

    def test_write_empty(self):
        s = _NullBinaryStream()
        assert s.write(b"") == 0
