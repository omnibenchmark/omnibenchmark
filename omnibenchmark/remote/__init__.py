"""
Remote storage integration for omnibenchmark.

Handles pushing and pulling benchmark artifacts to/from object storage so
results can be shared and versioned outside the local working tree.

Key modules:
- ``storage.py`` / ``RemoteStorage.py`` — storage abstraction and options
- ``S3Storage.py`` / ``S3config.py`` / ``S3versioning.py`` — S3-backed backend
  (also covers S3-compatible endpoints such as MinIO)
- ``files.py`` / ``tree.py`` / ``hash.py`` / ``sizeof.py`` — file listing,
  hashing, and size accounting for what gets transferred
- ``versioning.py`` — map a benchmark version to its expected output files
- ``archive.py`` / ``service.py`` — archive upload and high-level entry points

Layering: depends downward on :mod:`omnibenchmark.archive`,
:mod:`omnibenchmark.versioning`, and :mod:`omnibenchmark.config`. The remaining
``archive <-> remote`` cycle is the last one to untangle (a shared storage
abstraction would break it).
"""
