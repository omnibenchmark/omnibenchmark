"""
Exception classes for the versioning module.
"""


class VersioningError(Exception):
    """Base exception for all versioning-related errors."""

    pass


class VersionDowngradeError(VersioningError):
    """Raised when attempting to downgrade a version."""

    def __init__(self, current_version: str, proposed_version: str):
        self.current_version = current_version
        self.proposed_version = proposed_version
        super().__init__(
            f"Cannot downgrade from version {current_version} to {proposed_version}. "
            "Version downgrades are not allowed."
        )


class VersionLockError(VersioningError):
    """Raised when unable to acquire or release a version lock."""

    def __init__(self, lock_file: str, message: str = ""):
        self.lock_file = lock_file
        if message:
            super().__init__(f"Lock error for {lock_file}: {message}")
        else:
            super().__init__(f"Could not acquire lock for {lock_file}")


class VersionFormatError(VersioningError):
    """Raised when a version string has an invalid format."""

    def __init__(self, version_string: str, expected_format: str = "x.y.z or x.y"):
        self.version_string = version_string
        self.expected_format = expected_format
        super().__init__(
            f"Invalid version format: '{version_string}'. "
            f"Expected format: {expected_format}"
        )


class VersionAlreadyExistsError(VersioningError):
    """Raised when attempting to create a version that already exists."""

    def __init__(self, version: str):
        self.version = version
        super().__init__(f"Version {version} already exists")


class VersionNotFoundError(VersioningError):
    """Raised when a requested version cannot be found."""

    def __init__(self, version: str):
        self.version = version
        super().__init__(f"Version {version} not found")
