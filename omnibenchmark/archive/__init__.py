"""
Archive — bundle a benchmark's components into a portable archive.

Represents: assembling the parts of a benchmark (configuration, code, software
environments, results) into a single archive, for local-only use or upload.
When results come from object storage, the connected backend is *injected* by
the caller, so this package depends only on the storage interface.
Layer: service
Depends on: core, model, git, storage.
"""

from .archive import archive_benchmark
from .components import (
    prepare_archive_code,
    prepare_archive_software,
    prepare_archive_results,
    prepare_archive_config,
)

__all__ = [
    "archive_benchmark",
    "prepare_archive_code",
    "prepare_archive_software",
    "prepare_archive_results",
    "prepare_archive_config",
]
