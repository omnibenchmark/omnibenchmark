"""
Helper functions for building inputs_map in Snakemake workflows.

This module contains the logic for mapping template paths to actual file paths,
which is needed when Snakemake evaluates input files (especially for remote storage).
"""

from typing import Any, Dict


def build_inputs_map(
    inputs_map_template: Dict[str, Any], actual_input_files: list[str]
) -> Dict[str, Any]:
    """
    Build inputs_map by matching template paths to actual file paths.

    Snakemake passes inputs in the same order as defined in the rule, so we rely on
    positional matching. We iterate through actual files and assign them to keys
    in the order they appear in the template.

    For remote storage (S3), Snakemake downloads files to .snakemake/storage/ but
    preserves the order, so positional matching works correctly.

    Args:
        inputs_map_template: Template mapping from keys to expected paths
        actual_input_files: List of actual file paths from Snakemake (in same order)

    Returns:
        Mapping from keys to actual file paths

    Example:
        template = {
            "data.raw": "/tmp/D1.json",
            "methods.result": ["/tmp/M1.json", "/tmp/M2.json"]
        }
        actual = ["/storage/D1.json", "/storage/M1.json", "/storage/M2.json"]

        Returns: {
            "data.raw": "/storage/D1.json",
            "methods.result": ["/storage/M1.json", "/storage/M2.json"]
        }
    """
    inputs_map: Dict[str, Any] = {}
    file_iter = iter(actual_input_files)

    for key, value in inputs_map_template.items():
        if isinstance(value, str):
            # Single file - consume one from iterator
            inputs_map[key] = next(file_iter, None)
        elif isinstance(value, list):
            # Multiple files - consume len(value) files from iterator
            inputs_map[key] = [next(file_iter, None) for _ in value]

    return inputs_map
