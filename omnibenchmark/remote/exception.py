class RemoteStorageException(Exception):
    """
    General Exception for RemoteStorage.
    """

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


class RemoteStorageInvalidInputException(RemoteStorageException):
    """
    Exception raised for invalid inputs in RemoteStorage and descendants.
    """

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


class S3StorageException(RemoteStorageException):
    """
    General Exception for S3-compatible storage.
    """

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


class S3StorageConnectionException(S3StorageException):
    """
    Exception raised for connection issues with S3.
    """

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


class S3StorageBucketManipulationException(S3StorageException):
    """
    Exception raised for errors with bucket creation and general bucket manipulations.
    """

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


class S3StorageVersioningCorruptionException(S3StorageException):
    """
    Exception raised for errors with object versioning.
    """

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)
