"""
Storage abstraction — the contract for where benchmark artifacts live.

Represents: the *interface* (port) every storage backend implements, the value
objects that describe what gets tracked, and the pure rules for which files a
benchmark is expected to produce. It is the vocabulary of storage, with no
network access and no knowledge of any concrete backend.

Contains:
- ``base.RemoteStorage``    abstract backend interface (connect, download,
                            set_version, archive_version, ...)
- ``base.StorageOptions``   which directories/files count as tracked
- ``exceptions``            the storage exception hierarchy
- ``expected.get_expected_benchmark_output_files``  pure: benchmark + options
                            -> expected output file list (no I/O)

Layer: interface (port)
Depends on: ``core``.

It sits above :mod:`omnibenchmark.core` (it speaks in core domain objects) and
below every concrete or orchestrating layer.
Implemented by: :mod:`omnibenchmark.remote` (the concrete S3 backend + factory).
Consumed by: :mod:`omnibenchmark.remote`, :mod:`omnibenchmark.archive`.

Splitting this out of ``remote`` is what lets ``archive`` depend on the storage
*contract* without importing the storage *implementation* — breaking the old
``archive <-> remote`` import cycle.
"""

from .base import RemoteStorage, StorageOptions, is_valid_version
from .expected import get_expected_benchmark_output_files

__all__ = [
    "RemoteStorage",
    "StorageOptions",
    "is_valid_version",
    "get_expected_benchmark_output_files",
]
