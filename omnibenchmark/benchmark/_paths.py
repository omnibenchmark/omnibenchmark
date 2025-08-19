"""Path construction utilities for omnibenchmark workflows."""

import re
from pathlib import Path
from typing import List, Dict, Set, Optional


def normalize_dataset_path(path: str, prefix: Path) -> str:
    """Normalizes a benchmark output path by extracting and formatting the dataset identifier.

    Takes a path string containing a dataset identifier segment and formats it consistently.
    The path is expected to follow the structure: {prefix}/.../{dataset}/...

    Args:
        path: The output path string containing a dataset identifier
        prefix: The base directory Path object where outputs are stored

    Returns:
        A formatted path string with the dataset identifier properly substituted.
        Returns empty string if the path doesn't match the expected structure.

    Example:
        >>> normalize_dataset_path("out/stage1/dataset123/results", Path("out"))
        'out/stage1/dataset123/results'
    """
    pattern = rf"{prefix.as_posix()}/.+?/([^/]+)/.+?$"
    matched = re.match(pattern, path)
    if matched is None:
        # we probably should raise a ValueError here, but
        # this is better than an IndexError
        return ""
    dataset = matched[1]
    return path.format(dataset=dataset)


def construct_output_paths(
    prefix: Path,
    nodes: List,
    model,
    stage_outputs_cache: Optional[Dict[str, Dict]] = None,
) -> List[str]:
    """Construct output paths for a sequence of benchmark nodes.

    Args:
        prefix: Base directory path for outputs
        nodes: List of benchmark nodes representing execution path
        model: Benchmark model containing stage and module information
        stage_outputs_cache: Optional cache of stage outputs to avoid repeated lookups

    Returns:
        List of formatted output path strings
    """
    if nodes is None or len(nodes) == 0:
        return []

    head = nodes[0]
    tail = nodes[1:]

    # Use cache if provided, otherwise get stage outputs
    if stage_outputs_cache and head.stage_id in stage_outputs_cache:
        stage_outputs = stage_outputs_cache[head.stage_id].values()
    else:
        stage_outputs = model.get_stage_outputs(head.stage_id).values()

    current_path = f"{head.stage_id}/{head.module_id}"
    if any(["{params}" in o for o in stage_outputs]):
        current_path += f"/{head.param_id}"

    new_prefix = f"{prefix}/{current_path}"
    paths = [
        x.format(
            input=str(prefix),
            stage=head.stage_id,
            module=head.module_id,
            params=head.param_id,
            dataset="{dataset}",
        )
        for x in stage_outputs
    ]

    return paths + construct_output_paths(
        Path(new_prefix), tail, model, stage_outputs_cache
    )


def collect_output_paths(execution_paths: List[List], out_dir: Path, model) -> Set[str]:
    """Collect all output paths from execution paths.

    Args:
        execution_paths: List of execution paths, each containing a sequence of nodes
        out_dir: Output directory base path
        model: Benchmark model for stage output information

    Returns:
        Set of normalized output path strings
    """
    # Pre-cache stage outputs to avoid repeated lookups
    stage_outputs_cache = {}
    for stage_id, stage in model.get_stages().items():
        stage_outputs_cache[stage_id] = model.get_stage_outputs(stage)

    output_paths = []
    for path in execution_paths:
        path_outputs = construct_output_paths(
            prefix=out_dir,
            nodes=path,
            model=model,
            stage_outputs_cache=stage_outputs_cache,
        )
        for output in path_outputs:
            normalized_path = normalize_dataset_path(output, out_dir)
            if normalized_path:  # Only add non-empty paths
                output_paths.append(normalized_path)

    return set(output_paths)


def collect_path_exclusions(model) -> Dict[str, List[str]]:
    """Collect path exclusions from all modules in the benchmark.

    Args:
        model: Benchmark model containing stage and module information

    Returns:
        Dictionary mapping module IDs to their excluded module lists
    """
    path_exclusions = {}
    stages = model.get_stages()

    for stage_id in stages:
        stage = stages[stage_id]
        modules_in_stage = model.get_modules_by_stage(stage)

        for module_id in modules_in_stage:
            module_excludes = model.get_module_excludes(module_id)
            if module_excludes:
                path_exclusions[module_id] = module_excludes

    return path_exclusions
