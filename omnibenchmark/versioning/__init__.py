"""
Versioning Module for omnibenchmark.

This module consolidates all versioning logic for the omnibenchmark system into a single, cohesive
component. This design prevents versioning concerns from being scattered across
multiple modules and provides a unified foundation for all future versioning
requirements.

ARCHITECTURAL RATIONALE:
========================

The primary motivation for this centralized versioning module is architectural cleanliness.
Rather than having versioning logic distributed across storage layers, CLI commands, I/O
operations, and model classes, ALL versioning concerns are consolidated here. This provides:

1. **Single Source of Truth**: All version-related logic lives in one place
2. **Consistent Behavior**: All components use the same versioning rules and validation
3. **Future Extensibility**: New versioning features only need to be added here
4. **Easier Maintenance**: Version logic changes don't ripple across the entire codebase
5. **Clear Separation of Concerns**: Other modules focus on their primary responsibilities

MODULE'S APPROACH:
=========================

Instead of having things like strage clases doing their own version validation,
CLI commands scattering version comparison logic, or I/O operations with ad-hoc
version handling, or model classes mixing versioning with data representation...

This module allows to have:

- **ONE** module to rule them all: handles all version parsing, validation, and comparison
- **ONE** set of rules for version increments and constraints
- **ONE** interface for all version management operations
- **ONE** place to extend versioning capabilities in the future

ARCHITECTURAL LAYERS:
====================

1. **Core Version Logic** (version.py):
   - Version: Semantic version representation with comparison and increment operations
   - Utility functions for parsing and manipulating version strings
   - Base class on which other classes can build upon

2. **Version Management** (manager.py):
   - BenchmarkVersionManager: Core version lifecycle management with file-based
   locking. This allows for programmatic versioning, either for the CLI (instead
   of editing the YAML file directly) or from scripts, tests, etc. This can be
   done, for instance, to sync the git tag to a YAML version.
   - Concurrency control.
   - Hook system allowing integration without tight coupling to other modules
   (for instance, to auto-migrate fields when jumping across api versions).
   - In-memory state management with external history reconstruction

3. **Git Integration** (git.py):
   - GitAwareBenchmarkVersionManager: TENTATIVE implementation extending base manager
   - **STATUS**: Example implementation pending proper design specification. We need
     to define the API and how benchmark versioning interacts with remote
     repositories (e.g. S3 archiving)
   - Again: This is exploratory code demonstrating potential git integration
     patterns. It can be extended or trown away.
   - Shows how the centralized architecture can accommodate different versioning backends

4. **Exception Hierarchy** (exceptions.py):
   - Unified exception types for all versioning error scenarios across the system
   - Prevents each module from defining its own version-related exceptions

DESIGN PRINCIPLES:
=================

- **Centralization**: All versioning logic consolidated in one place
- **Extensibility**: New versioning strategies can be added without affecting existing code
- **Consistency**: Same versioning behavior across all system components
- **Separation**: Other modules delegate versioning concerns to this module

MAIN CLASSES:
============

Version:
    Core semantic version representation (x.y.z or x.y format).
    Handles all version parsing, comparison, and increment operations.
    Used by all other components requiring version operations.

BenchmarkVersionManager:
    Version manager providing:
    - File-based locking for concurrency safety
    - Version validation and lifecycle management
    - Hook system for extensible integration
    - Monotonic versioning with downgrade prevention

GitAwareBenchmarkVersionManager:
    **TENTATIVE IMPLEMENTATION** - Example of extending base manager.
    Shows potential patterns for:
    - Version history reconstruction from git commits
    - Automatic benchmark file updates with git tracking
    - Git repository state integration

    **IMPORTANT**: This is exploratory code demonstrating how the module
    can accommodate git-based versioning. The implementation is
    incomplete and serves as a design example pending proper specification.

SYSTEM INTEGRATION:
==================

This centralized versioning module serves as the foundation for ALL versioning
operations across the omnibenchmark system:

1. **Storage Layer** (omnibenchmark.io.*):
   - Delegates ALL version validation to this module
   - Uses version managers for storage synchronization
   - No storage-specific versioning logic needed

2. **CLI Commands** (omnibenchmark.cli.*):
   - All version operations route through this module
   - Consistent version behavior across all commands
   - No CLI-specific version handling needed

3. **I/O Operations**:
   - Version tagging and correlation handled here
   - File versioning operations centralized
   - No scattered I/O version logic

4. **Model Layer**:
   - Models focus on data representation only
   - Version validation delegated to this module
   - Clean separation between data and versioning logic

FUTURE EXTENSIBILITY:
=====================

A centralized module makes future versioning enhancements straightforward:

- **New Version Schemes**: Add new Version subclasses without changing other modules
- **Different Backends**: Extend managers for database, cloud, or other version stores
- **Advanced Features**: Add branching, tagging, or complex workflows in one place
- **Integration Points**: New systems integrate through the established interfaces

All future versioning requirements can be accommodated by extending this module
rather than adding versioning logic throughout the codebase.
"""

from .manager import BenchmarkVersionManager
from .git import GitAwareBenchmarkVersionManager
from .exceptions import (
    VersioningError,
    VersionDowngradeError,
    VersionLockError,
    VersionFormatError,
    VersionAlreadyExistsError,
    VersionNotFoundError,
)
from .version import Version, parse_version, increment_version

__all__ = [
    # Main production-ready classes
    "BenchmarkVersionManager",
    # Tentative/example implementations
    "GitAwareBenchmarkVersionManager",  # TENTATIVE - example pending design spec
    # Core version utilities
    "Version",
    "parse_version",
    "increment_version",
    # Centralized exception hierarchy
    "VersioningError",
    "VersionDowngradeError",
    "VersionLockError",
    "VersionFormatError",
    "VersionAlreadyExistsError",
    "VersionNotFoundError",
]
