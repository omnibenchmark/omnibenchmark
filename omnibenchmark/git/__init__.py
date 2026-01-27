"""
Git operations module for omnibenchmark.

This module provides functionality for cloning and managing git repositories
used in benchmark workflows, separated from remote storage concerns.
"""

from .clone import clone_git_repo, clone_module

__all__ = ["clone_git_repo", "clone_module"]
