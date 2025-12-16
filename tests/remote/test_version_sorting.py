"""
Unit test to prevent regression of version sorting bug.

This test guards against the TypeError that occurred when attempting to use
packaging.version.Version as a key function instead of sorting Version objects directly.
"""

import pytest
from packaging.version import Version


@pytest.mark.short
def test_version_sorting_regression():
    """
    Test that Version objects can be sorted directly without using Version as key function.

    This guards against the regression where code attempted:
        versions.sort(key=Version)  # WRONG - TypeError

    Instead of the correct:
        versions.sort()  # CORRECT - Version objects are comparable
    """
    # Create a list of Version objects (similar to what remote version list does)
    version_strings = ["2.0", "1.0", "1.5", "0.9", "10.0"]
    versions = [Version(v) for v in version_strings]

    # This should work - Version objects are directly comparable
    versions_copy = versions.copy()
    versions_copy.sort()

    # Verify they're sorted correctly
    expected_order = [
        Version("0.9"),
        Version("1.0"),
        Version("1.5"),
        Version("2.0"),
        Version("10.0"),
    ]
    assert versions_copy == expected_order

    # This would cause TypeError if attempted (what the bug was doing):
    # versions.sort(key=Version)  # Don't do this!

    # Verify the bug scenario would indeed fail
    # The original bug was: versions.sort(key=Version) where versions is a list of Version objects
    # This fails because Version class itself can't be used as a key function
    with pytest.raises(TypeError):
        # This simulates the original buggy code
        version_objects = [Version("1.0"), Version("2.0")]
        version_objects.sort(
            key=Version
        )  # This should fail - Version class used as key function


@pytest.mark.short
def test_version_objects_are_comparable():
    """Test that Version objects implement comparison operators correctly."""
    v1 = Version("1.0")
    v2 = Version("2.0")
    v3 = Version("1.0")

    # Test comparison operators
    assert v1 < v2
    assert v2 > v1
    assert v1 == v3
    assert v1 <= v2
    assert v2 >= v1
    assert v1 != v2

    # Test that they work in sorted()
    versions = [v2, v1, v3]
    sorted_versions = sorted(versions)
    assert sorted_versions == [v1, v3, v2]


@pytest.mark.short
def test_version_list_sorting_pattern():
    """
    Test the correct pattern for sorting version lists as used in remote commands.

    This demonstrates the correct way to sort a list of Version objects
    that would be returned by remote version list functionality.
    """

    # Simulate what remote version list might return
    class MockVersionList:
        def __init__(self, version_strings):
            self.versions = [Version(v) for v in version_strings]

        def sort_versions(self):
            # This is the CORRECT way (what the fix implemented)
            self.versions.sort()
            return self.versions

    # Test with realistic version numbers
    version_strings = ["1.0", "2.1", "1.5", "2.0", "1.0-beta"]
    mock_list = MockVersionList(version_strings)

    sorted_versions = mock_list.sort_versions()

    # Verify they're in correct order
    expected = [
        Version("1.0-beta"),
        Version("1.0"),
        Version("1.5"),
        Version("2.0"),
        Version("2.1"),
    ]
    assert sorted_versions == expected

    # Verify the original list was modified (sort in place)
    assert mock_list.versions == expected
