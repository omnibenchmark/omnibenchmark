"""Path construction utilities for omnibenchmark workflows."""

import hashlib
import re
from pathlib import PurePosixPath, Path
from typing import List, Dict, Set, Optional

# Per-component limit on most modern filesystems (ext4, xfs, apfs, ntfs).
MAX_FILENAME_LEN = 255

# Cap on what we treat as a "compound extension" before falling back to the
# last suffix only. Covers .tar.gz, .tar.bz2, .tar.zst, .nii.gz, etc.
_MAX_COMPOUND_SUFFIX_LEN = 16
_EXT_RE = re.compile(r"\.[A-Za-z0-9]+")


def _split_extension(name: str) -> tuple[str, str]:
    """Split *name* into (stem, suffix), preserving compound extensions.

    Only treats trailing dot-segments matching ``.[A-Za-z0-9]+`` as extension,
    and only when the joined suffix is at most ``_MAX_COMPOUND_SUFFIX_LEN``
    characters — so ``foo.tar.gz`` -> ("foo", ".tar.gz") but a pathological
    ``my.dotted.name.txt`` -> ("my.dotted.name", ".txt") rather than eating
    the whole stem.
    """
    suffixes = PurePosixPath(name).suffixes
    if not suffixes:
        return name, ""

    compound = "".join(suffixes)
    if len(compound) <= _MAX_COMPOUND_SUFFIX_LEN and all(
        _EXT_RE.fullmatch(s) for s in suffixes
    ):
        return name[: -len(compound)], compound

    last = suffixes[-1]
    if _EXT_RE.fullmatch(last) and len(last) <= _MAX_COMPOUND_SUFFIX_LEN:
        return name[: -len(last)], last

    return name, ""


def truncate_filename(name: str, limit: int = MAX_FILENAME_LEN) -> str:
    """Truncate a single filename component to fit *limit* bytes.

    Preserves the file extension (including compound forms like ``.tar.gz``)
    and appends an 8-char hash of the original stem for uniqueness, so two
    long names that share a prefix still produce distinct truncated names.

    Operates on a single path component — callers must split parent dirs off
    before calling.
    """
    if len(name) <= limit:
        return name

    stem, suffix = _split_extension(name)
    digest = hashlib.sha256(stem.encode("utf-8")).hexdigest()[:8]
    sep = "_"
    keep = limit - len(suffix) - len(sep) - len(digest)
    if keep < 1:
        # Limit too small to fit hash + suffix; degrade gracefully.
        return (digest + suffix)[:limit]
    return stem[:keep] + sep + digest + suffix


def sanitize_rule_name(node_id: str) -> str:
    """Sanitize a node ID into a valid Snakemake rule name (Python identifier).

    Shared by the Snakefile generator (which names each rule's ``log:`` file)
    and the status reporter (which locates those log files), so the two stay in
    lockstep.
    """
    name = node_id.replace("-", "_").replace(".", "_")
    if name and not name[0].isalpha():
        name = "rule_" + name
    return name


_UNSAFE_CHARS = ["/", "\\", ":", "*", "?", '"', "<", ">", "|", " "]


def make_human_name(params, limit: int = MAX_FILENAME_LEN) -> str:
    """Create a human-readable symlink name from a Params object.

    Joins key-value pairs with underscores and strips filesystem-unsafe
    characters.  Falls back to appending the params' short hash if the result
    exceeds *limit* (the typical filesystem per-component limit).
    """
    parts = [f"{k}-{v}" for k, v in params.items()]
    name = "_".join(parts)
    for char in _UNSAFE_CHARS:
        name = name.replace(char, "_")
    if len(name) > limit:
        digest = params.hash_short()
        keep = limit - len(digest) - 1
        name = name[:keep] + "_" + digest
    return name


def truncate_path_filename(path: str, limit: int = MAX_FILENAME_LEN) -> str:
    """Apply :func:`truncate_filename` to the leaf component of *path* only."""
    parent, sep, basename = path.rpartition("/")
    truncated = truncate_filename(basename, limit=limit)
    if truncated == basename:
        return path
    return f"{parent}{sep}{truncated}" if sep else truncated


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


def is_lineage_excluded(
    module_ids: Set[str], path_exclusions: Dict[str, List[str]]
) -> bool:
    """Whether a lineage of module IDs violates any declared exclusion.

    ``module_ids`` is the full set of module IDs along a path (e.g. an execution
    path or a candidate node's lineage). A lineage is excluded when, for any
    module that declares exclusions, both the declaring module and one of its
    excluded modules are present — i.e. the degenerate combination co-occurs on
    the same path, regardless of how many stages separate them.

    This is the single source of truth shared by execution-path pruning
    (``exclude_paths``) and Snakefile generation (``cli/run.py``).
    """
    for module_id, excluded in path_exclusions.items():
        if module_id in module_ids and any(e in module_ids for e in excluded):
            return True
    return False
