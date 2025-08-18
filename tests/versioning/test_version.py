"""
Tests for the Version class and version utilities.

All tests in this file are marked as 'short' since they don't require
external dependencies, containers, or network I/O.
"""

import pytest
from omnibenchmark.versioning.version import (
    Version,
    parse_version,
    increment_version,
    compare_versions,
)


@pytest.mark.short
class TestVersion:
    """Test the Version class."""

    def test_version_creation_valid(self):
        """Test creating versions with valid formats."""
        # Test x.y.z format
        v1 = Version("1.2.3")
        assert v1.major == 1
        assert v1.minor == 2
        assert v1.patch == 3
        assert str(v1) == "1.2.3"

        # Test x.y format
        v2 = Version("2.1")
        assert v2.major == 2
        assert v2.minor == 1
        assert v2.patch is None
        assert str(v2) == "2.1"

        # Test with zeros
        v3 = Version("0.0.0")
        assert v3.major == 0
        assert v3.minor == 0
        assert v3.patch == 0
        assert str(v3) == "0.0.0"

    def test_version_creation_invalid(self):
        """Test that invalid version formats raise ValueError."""
        invalid_versions = [
            "1",  # Missing minor
            "a.b.c",  # Non-numeric
            "1.2.3.4",  # Too many parts
            "1.2.a",  # Non-numeric patch
            "1.-2",  # Negative number
            "",  # Empty string
            "v1.2.3",  # With prefix
        ]

        for invalid in invalid_versions:
            with pytest.raises(ValueError, match="Invalid version format"):
                Version(invalid)

    def test_version_equality(self):
        """Test version equality comparison."""
        v1 = Version("1.2.3")
        v2 = Version("1.2.3")
        v3 = Version("1.2")
        v4 = Version("1.2.0")

        assert v1 == v2
        assert v1 != v3
        assert v3 != v4  # 1.2 != 1.2.0 (None != 0)

        # Test with non-Version object
        assert v1 != "1.2.3"
        assert v1 != 123
        assert v1 is not None

    def test_version_comparison(self):
        """Test version comparison operators."""
        v1 = Version("1.0.0")
        v2 = Version("1.0.1")
        v3 = Version("1.1.0")
        v4 = Version("2.0.0")

        # Less than
        assert v1 < v2
        assert v2 < v3
        assert v3 < v4

        # Greater than
        assert v4 > v3
        assert v3 > v2
        assert v2 > v1

        # Less than or equal
        assert v1 <= v1
        assert v1 <= v2

        # Greater than or equal
        assert v2 >= v2
        assert v2 >= v1

    def test_version_comparison_with_missing_patch(self):
        """Test version comparison when patch version is missing."""
        v1 = Version("1.2")
        v2 = Version("1.2.0")
        v3 = Version("1.2.1")

        # 1.2 (None) is treated as 0 for comparison
        assert v1 < v3
        assert v2 < v3

    def test_version_hash(self):
        """Test that versions can be used in sets and as dict keys."""
        v1 = Version("1.2.3")
        v2 = Version("1.2.3")
        v3 = Version("1.2.4")

        # Same versions should have same hash
        assert hash(v1) == hash(v2)

        # Can be used in sets
        version_set = {v1, v2, v3}
        assert len(version_set) == 2  # v1 and v2 are the same

        # Can be used as dict keys
        version_dict = {v1: "first", v3: "second"}
        assert version_dict[v2] == "first"  # v2 equals v1

    def test_version_increment_minor(self):
        """Test incrementing minor version."""
        v1 = Version("1.2.3")
        v2 = v1.increment_minor()
        assert str(v2) == "1.3"

        v3 = Version("2.0")
        v4 = v3.increment_minor()
        assert str(v4) == "2.1"

    def test_version_increment_major(self):
        """Test incrementing major version."""
        v1 = Version("1.2.3")
        v2 = v1.increment_major()
        assert str(v2) == "2.0"

        v3 = Version("0.9")
        v4 = v3.increment_major()
        assert str(v4) == "1.0"

    def test_version_increment_patch(self):
        """Test incrementing patch version."""
        v1 = Version("1.2.3")
        v2 = v1.increment_patch()
        assert str(v2) == "1.2.4"

        # When patch is None, start at 1
        v3 = Version("1.2")
        v4 = v3.increment_patch()
        assert str(v4) == "1.2.1"

    def test_version_repr(self):
        """Test string representation of Version."""
        v1 = Version("1.2.3")
        assert repr(v1) == "Version('1.2.3')"

        v2 = Version("2.1")
        assert repr(v2) == "Version('2.1')"


@pytest.mark.short
class TestVersionUtilities:
    """Test version utility functions."""

    def test_parse_version(self):
        """Test parse_version function."""
        v1 = parse_version("1.2.3")
        assert isinstance(v1, Version)
        assert str(v1) == "1.2.3"

        v2 = parse_version("2.1")
        assert isinstance(v2, Version)
        assert str(v2) == "2.1"

        with pytest.raises(ValueError):
            parse_version("invalid")

    def test_increment_version(self):
        """Test increment_version function."""
        # Test minor increment (default)
        assert increment_version("1.2.3") == "1.3"
        assert increment_version("1.2") == "1.3"
        assert increment_version("1.2.3", "minor") == "1.3"

        # Test major increment
        assert increment_version("1.2.3", "major") == "2.0"
        assert increment_version("0.9", "major") == "1.0"

        # Test patch increment
        assert increment_version("1.2.3", "patch") == "1.2.4"
        assert increment_version("1.2", "patch") == "1.2.1"

        # Test invalid component
        with pytest.raises(ValueError, match="Unknown version component"):
            increment_version("1.2.3", "invalid")

        # Test invalid version
        with pytest.raises(ValueError):
            increment_version("invalid", "minor")

    def test_compare_versions(self):
        """Test compare_versions function."""
        # Equal versions
        assert compare_versions("1.2.3", "1.2.3") == 0
        assert compare_versions("1.2", "1.2") == 0

        # First less than second
        assert compare_versions("1.2.3", "1.2.4") == -1
        assert compare_versions("1.2", "2.0") == -1
        assert compare_versions("0.9", "1.0") == -1

        # First greater than second
        assert compare_versions("1.2.4", "1.2.3") == 1
        assert compare_versions("2.0", "1.2") == 1
        assert compare_versions("1.0", "0.9") == 1

        # Test with invalid versions
        with pytest.raises(ValueError):
            compare_versions("invalid", "1.2.3")

        with pytest.raises(ValueError):
            compare_versions("1.2.3", "invalid")


@pytest.mark.short
class TestVersionEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_large_version_numbers(self):
        """Test with large version numbers."""
        v1 = Version("999.999.999")
        assert v1.major == 999
        assert v1.minor == 999
        assert v1.patch == 999

        v2 = v1.increment_patch()
        assert str(v2) == "999.999.1000"

    def test_version_sorting(self):
        """Test that versions sort correctly."""
        versions = [
            Version("2.0.0"),
            Version("1.0.0"),
            Version("1.2.0"),
            Version("1.1.0"),
            Version("1.1.1"),
            Version("1.1"),
        ]

        sorted_versions = sorted(versions)

        # Our implementation treats None patch as 0 for comparison
        # but they are not equal (1.1 != 1.1.0 in __eq__)
        # Since they compare as equal in __lt__, stable sort preserves input order
        expected = [
            "1.0.0",
            "1.1.0",  # Was at index 3 in input
            "1.1",  # Was at index 5 in input (comes after 1.1.0 due to stable sort)
            "1.1.1",
            "1.2.0",
            "2.0.0",
        ]

        assert [str(v) for v in sorted_versions] == expected

    def test_version_whitespace_handling(self):
        """Test that whitespace is properly handled."""
        v1 = Version("  1.2.3  ")
        assert str(v1) == "1.2.3"

        v2 = Version("\t2.1\n")
        assert str(v2) == "2.1"

    def test_version_string_conversion(self):
        """Test string conversion edge cases."""
        # Float with trailing zeros
        v1 = Version("1.0")
        assert str(v1) == "1.0"

        # Integer string
        v2 = Version("2.0")
        assert str(v2) == "2.0"
