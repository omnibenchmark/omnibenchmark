"""CLI commands for collecting metrics from a benchmark output folder.

The ``collect`` group gathers post-hoc artefacts that are scattered across the
output tree into a single combined table. Right now it knows how to gather the
Snakemake ``performance.txt`` benchmark files emitted for every executed node;
the resulting ``performances.tsv`` is what ``ob dashboard`` consumes.

Collection is a pure filesystem walk -- it does not invoke Snakemake and can be
run against any existing output directory.
"""

import csv
import json
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import click

from omnibenchmark.constants import INTERNAL_OUT_DIRS
from omnibenchmark.logging import logger
from .debug import add_debug_option

# Internal state dirs, plus `.metadata` and the archived `old` runs: never a
# source of performance files, expensive (and noisy) to walk.
_SKIP_DIRS = INTERNAL_OUT_DIRS | {".metadata", "old"}

# A node's benchmark file is `performance.txt` (>= v0.5.0) or, in older layouts,
# `{dataset}_performance.txt`. Match both so the collector works across versions.
_PERF_SUFFIX = "_performance.txt"
_PERF_NAME = "performance.txt"

# Column emitted by Snakemake that is redundant with `s` and not numeric.
_DROP_COLS = {"h:m:s"}


def _is_performance_file(name: str) -> bool:
    return name == _PERF_NAME or name.endswith(_PERF_SUFFIX)


def _find_performance_files(out_dir: Path) -> List[Path]:
    """Walk ``out_dir`` and return every performance file, skipping noise dirs."""
    found: List[Path] = []
    for root, dirs, files in os.walk(out_dir):
        # Prune skip dirs in place so os.walk does not descend into them.
        dirs[:] = [d for d in dirs if d not in _SKIP_DIRS]
        for name in files:
            if _is_performance_file(name):
                found.append(Path(root) / name)
    return sorted(found)


def _read_performance_row(file_path: Path) -> Optional[Dict[str, Any]]:
    """Read the single data row of a Snakemake benchmark TSV, coercing numbers."""
    with open(file_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for raw in reader:
            row: Dict[str, Any] = {}
            for key, value in raw.items():
                if key in _DROP_COLS:
                    continue
                if value is None or value in ("", "NA", "none", "null", "nan"):
                    row[key] = 0
                else:
                    try:
                        row[key] = float(value)
                    except (ValueError, TypeError):
                        row[key] = value
            return row  # benchmark files hold exactly one data row
    return None


def _split_gather_group(rel_parts: List[str]) -> tuple:
    """Pull the group key out of a grouped gather's path, if there is one.

    A grouped gather roots at ``stage/<group>/module/<param-dir>`` (design 010
    §3.3), four segments where every other node contributes three. So a single
    leftover component makes the path self-describing: it is the group key, in
    second position, and removing it restores the 3-aligned chain.
    """
    parts = list(rel_parts)
    if len(parts) > 1 and len(parts) % 3 == 1:
        return parts.pop(1), parts
    return None, parts


def _triples(rel_parts: List[str]) -> List[tuple]:
    """Group path components into (stage, module, param_dir) triples.

    The output layout is a chain of ``stage/module/<param-dir>`` directories,
    where ``<param-dir>`` is ``.default`` (or ``default`` in old layouts) for
    parameter-free nodes and ``.<hash>`` otherwise. Any gather group key must
    already be removed (see ``_split_gather_group``).
    """
    return list(zip(*(iter(rel_parts),) * 3))


def _is_default_param_dir(param_dir: str) -> bool:
    return param_dir.lstrip(".") == "default"


def _collect_params(
    out_dir: Path, triples: List[tuple], group: Optional[str] = None
) -> Dict[str, Any]:
    """Read every ``parameters.json`` along the lineage into ``{stage: params}``."""
    params: Dict[str, Any] = {}
    cursor = out_dir
    for index, (stage, module, param_dir) in enumerate(triples):
        cursor = cursor / stage
        if index == 0 and group is not None:
            cursor = cursor / group
        cursor = cursor / module / param_dir
        if _is_default_param_dir(param_dir):
            continue
        param_file = cursor / "parameters.json"
        if param_file.exists():
            try:
                params[stage] = json.loads(param_file.read_text())
            except (json.JSONDecodeError, OSError):
                logger.warning(f"Could not parse parameters at {param_file}")
    return params


def _build_record(out_dir: Path, perf_file: Path) -> Optional[Dict[str, Any]]:
    """Turn one performance file into an enriched row, or None if unreadable."""
    row = _read_performance_row(perf_file)
    if row is None:
        logger.warning(f"No performance data in {perf_file}; skipping.")
        return None

    rel_parts = perf_file.relative_to(out_dir).parts[:-1]  # drop filename
    group, chain = _split_gather_group(list(rel_parts))
    triples = _triples(chain)
    if not triples:
        logger.warning(f"Unexpected path layout for {perf_file}; skipping.")
        return None

    leaf_stage, leaf_module, leaf_param_dir = triples[-1]
    params = _collect_params(out_dir, triples, group)

    # Metadata columns; names match what `ob dashboard` excludes from metrics.
    row["stage"] = leaf_stage
    row["module"] = leaf_module
    # Below a gather the chain is cut, so the group key is the only thing left
    # identifying the data the row came from.
    row["dataset"] = group if group is not None else triples[0][1]
    row["param_hash"] = (
        "" if _is_default_param_dir(leaf_param_dir) else leaf_param_dir.lstrip(".")
    )
    row["params"] = json.dumps(params) if params else ""
    row["lineage"] = "/".join(f"{s}/{m}" for s, m, _ in triples)
    row["path"] = str(perf_file.relative_to(out_dir))
    return row


def _write_tsv(output_file: Path, rows: List[Dict[str, Any]]) -> None:
    # Stable column order: metrics first (from the first row), metadata appended
    # in a fixed order so the file is diff-friendly across runs.
    meta_cols = [
        "stage",
        "module",
        "dataset",
        "param_hash",
        "params",
        "lineage",
        "path",
    ]
    metric_cols = [c for c in rows[0].keys() if c not in meta_cols]
    fieldnames = metric_cols + meta_cols
    with open(output_file, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


@click.group(name="collect")
def collect():
    """Gather scattered benchmark artefacts into combined tables."""


@add_debug_option
@collect.command(name="performance")
@click.option(
    "--out-dir",
    "-o",
    type=click.Path(),
    default="out",
    help="Output directory containing benchmark results (default: out).",
)
@click.pass_context
def performance(ctx, out_dir):
    """Gather all performance.txt files into a combined performances.tsv.

    Walks the output directory, parses every Snakemake benchmark file, and
    enriches each row with its stage, module, parameters and lineage. The
    resulting out/performances.tsv is consumed by `ob dashboard`.
    """
    out_path = Path(out_dir)
    if not out_path.exists():
        logger.error(f"Output directory does not exist: {out_dir}")
        sys.exit(1)

    perf_files = _find_performance_files(out_path)
    if not perf_files:
        logger.error(f"No performance files found under {out_dir}.")
        sys.exit(1)

    rows: List[Dict[str, Any]] = []
    for perf_file in perf_files:
        record = _build_record(out_path, perf_file)
        if record is not None:
            rows.append(record)

    if not rows:
        logger.error("No readable performance data collected.")
        sys.exit(1)

    output_file = out_path / "performances.tsv"
    _write_tsv(output_file, rows)
    logger.info(f"Collected {len(rows)} performance records to {output_file}.")
    sys.exit(0)
