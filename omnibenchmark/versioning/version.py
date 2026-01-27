"""
Version utility module for version string operations.

This module provides utilities for working with semantic versions,
using the standard packaging.version library for robust version handling.
"""

from typing import Optional
import re
from packaging.version import Version as PackagingVersion, InvalidVersion


class Version:
    """
    A semantic version representation using packaging.version.

    This class wraps packaging.version.Version to provide increment operations
    and maintain compatibility with the existing API.
    Version format: x.y.z or x.y (where x, y, z are non-negative integers)
    """

    def __init__(self, version_string: str):
        """
        Initialize a Version from a string.

        Args:
            version_string: Version string in format "x.y.z" or "x.y"

        Raises:
            ValueError: If version string is invalid
        """
        # Handle numeric inputs (float/int from YAML)
        if isinstance(version_string, (int, float)):
            version_string = str(version_string)

        self._original_string = str(version_string).strip()

        # Validate strict semantic version format first
        pattern = r"^(\d+)\.(\d+)(?:\.(\d+))?$"
        if not re.match(pattern, self._original_string):
            raise ValueError(
                f"Invalid version format: '{self._original_string}'. Expected x.y.z or x.y"
            )

        try:
            self._version = PackagingVersion(self._original_string)

            # Additional validation for pre/post/dev releases
            if (
                self._version.is_prerelease
                or self._version.is_postrelease
                or self._version.is_devrelease
            ):
                raise ValueError(
                    f"Only release versions are supported: '{self._original_string}'"
                )

        except InvalidVersion as e:
            raise ValueError(
                f"Invalid version format: '{self._original_string}'. Expected x.y.z or x.y"
            ) from e

    @property
    def major(self) -> int:
        """Major version component."""
        return self._version.major

    @property
    def minor(self) -> int:
        """Minor version component."""
        return self._version.minor

    @property
    def patch(self) -> Optional[int]:
        """Patch version component (None if not specified in original string)."""
        # Check if original string had patch component
        if self._original_string.count(".") >= 2:
            return self._version.micro
        return None

    def __str__(self) -> str:
        """Return the string representation of the version."""
        if self.patch is None:
            return f"{self.major}.{self.minor}"
        return f"{self.major}.{self.minor}.{self.patch}"

    def __repr__(self) -> str:
        """Return the debug representation of the version."""
        return f"Version('{str(self)}')"

    def __eq__(self, other) -> bool:
        """Check if two versions are equal."""
        if not isinstance(other, Version):
            return False
        # Use original format comparison to preserve 1.2 != 1.2.0 behavior
        return (self.major, self.minor, self.patch) == (
            other.major,
            other.minor,
            other.patch,
        )

    def __lt__(self, other) -> bool:
        """Check if this version is less than another."""
        if not isinstance(other, Version):
            return NotImplemented
        return self._version < other._version

    def __le__(self, other) -> bool:
        """Check if this version is less than or equal to another."""
        if not isinstance(other, Version):
            return NotImplemented
        return self._version <= other._version

    def __gt__(self, other) -> bool:
        """Check if this version is greater than another."""
        if not isinstance(other, Version):
            return NotImplemented
        return self._version > other._version

    def __ge__(self, other) -> bool:
        """Check if this version is greater than or equal to another."""
        if not isinstance(other, Version):
            return NotImplemented
        return self._version >= other._version

    def __hash__(self) -> int:
        """Return hash of the version for use in sets/dicts."""
        return hash(self._version)

    def increment_minor(self) -> "Version":
        """Return a new Version with incremented minor version."""
        return Version(f"{self.major}.{self.minor + 1}")

    def increment_major(self) -> "Version":
        """Return a new Version with incremented major version."""
        return Version(f"{self.major + 1}.0")

    def increment_patch(self) -> "Version":
        """Return a new Version with incremented patch version."""
        patch = self.patch if self.patch is not None else 0
        return Version(f"{self.major}.{self.minor}.{patch + 1}")


def parse_version(version_string: str) -> Version:
    """
    Parse a version string into a Version object.

    Args:
        version_string: Version string to parse

    Returns:
        Version object

    Raises:
        ValueError: If version string is invalid
    """
    return Version(version_string)


def increment_version(version: str, component: str = "minor") -> str:
    """
    Increment a version string.

    Args:
        version: Current version string
        component: Which component to increment ("major", "minor", or "patch")

    Returns:
        Incremented version string

    Raises:
        ValueError: If version string is invalid or component is unknown
    """
    v = Version(version)

    if component == "major":
        new_v = v.increment_major()
    elif component == "minor":
        new_v = v.increment_minor()
    elif component == "patch":
        new_v = v.increment_patch()
    else:
        raise ValueError(f"Unknown version component: {component}")

    return str(new_v)


def compare_versions(version1: str, version2: str) -> int:
    """
    Compare two version strings.

    Args:
        version1: First version string
        version2: Second version string

    Returns:
        -1 if version1 < version2
         0 if version1 == version2
         1 if version1 > version2

    Raises:
        ValueError: If either version string is invalid
    """
    v1 = Version(version1)
    v2 = Version(version2)

    if v1 < v2:
        return -1
    elif v1 > v2:
        return 1
    else:
        return 0
