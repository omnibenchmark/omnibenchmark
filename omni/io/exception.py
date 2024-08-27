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


class MinIOStorageException(RemoteStorageException):
    """
    General Exception for MinIOStorage.
    """

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


class MinIOStorageConnectionException(MinIOStorageException):
    """
    Exception raised for connection issues with S3.
    """

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


class MinIOStorageBucketManipulationException(MinIOStorageException):
    """
    Exception raised for errors with bucket creation and general bucket manipulations.
    """

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


class MinIOStorageVersioningCorruptionException(MinIOStorageException):
    """
    Exception raised for errors with object versioning.
    """

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)
