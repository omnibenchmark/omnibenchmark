"""
Remote storage backend — the concrete implementation of the storage interface.

Represents: actually moving bytes to/from object storage so results can be
shared and versioned outside the local working tree. This is the *wire*; the
*contract* it implements lives in :mod:`omnibenchmark.storage`.
Layer: backend (implementation)
Depends on: storage, archive, versioning, core, config.

Key modules:
- ``S3Storage.py`` / ``S3config.py`` / ``S3versioning.py`` — S3-backed backend
  (also covers S3-compatible endpoints such as MinIO)
- ``storage.py`` — the StorageFactory / get_storage backend selection, plus
  credential and Snakemake-arg helpers
- ``files.py`` / ``tree.py`` / ``hash.py`` / ``sizeof.py`` — file listing,
  hashing, and size accounting for what gets transferred
- ``versioning.py`` — object tagging/filtering for benchmark versions
- ``archive.py`` / ``service.py`` — archive upload and high-level entry points
"""
