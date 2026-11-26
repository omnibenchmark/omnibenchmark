"""Which benchmark an output dir belongs to, and which part of it each run built.

`ob run` appends one line per invocation to `out/.metadata/runs.jsonl` (007 §4.6):
the benchmark identity (`id`, `version`, `summary_hash`) and the selection. The
readers (`ob validate outputs`, `ob describe status`) compare the plan they were
given against the recorded identity, not against a hash recomputed from the
plan copy (a newer ob may hash the same YAML differently). Only an out/ written
before runs.jsonl existed falls back to the copy.
"""

import hashlib
import json
import shutil
from pathlib import Path
from typing import Iterable, List, Optional, Tuple

import yaml

from omnibenchmark.logging import logger

METADATA_DIR = ".metadata"
RUNLOG = "runs.jsonl"
PLAN_COPY = "benchmark.yaml"


def _sha256(obj) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True).encode()).hexdigest()


def benchmark_record(model) -> dict:
    return {
        "id": model.id,
        "version": model.get_version(),
        "summary_hash": model.summary_hash(),
        "api_version": model.api_version.value,
    }


def selection_of(
    module: Optional[str], until: Optional[str], filter_blob: Optional[dict]
) -> Optional[dict]:
    """None for a full run, else what restricted it. The three are exclusive."""
    if module:
        return {"module": module}
    if until:
        return {"until": until}
    if filter_blob is not None:
        return {
            "filter": {
                "sha256": _sha256(filter_blob),
                "parent": (filter_blob.get("parent") or {}).get("sha256"),
            }
        }
    return None


def write_run_record(
    out_dir: Path,
    manifest: dict,
    plan_path: Path,
    node_ids: Iterable[str],
    selection: Optional[dict],
    filter_blob: Optional[dict] = None,
    capabilities: Iterable[str] = (),
    trace: Optional[str] = None,
) -> dict:
    """Append this run to runs.jsonl, next to content-addressed plan, filter and node list.

    *manifest* is what `write_run_manifest` returned: run_id, timestamp,
    ob_version and the benchmark identity are taken from it, so the two records
    cannot disagree about the same run. With telemetry on, run_id is the trace's
    traceId and *trace* says where the spans went (a file path, or "stdout").
    """
    md = Path(out_dir) / METADATA_DIR
    md.mkdir(parents=True, exist_ok=True)
    bench = manifest["benchmark"]

    plan_copy = md / f"benchmark-{bench['summary_hash'][:8]}.yaml"
    if not plan_copy.exists():
        shutil.copy(plan_path, plan_copy)
    if selection and "filter" in selection and filter_blob is not None:
        sha = selection["filter"]["sha256"]
        (md / f"filter-{sha[:8]}.json").write_text(
            json.dumps(filter_blob, indent=2, sort_keys=True) + "\n"
        )
    nodes = sorted(node_ids)
    nodes_file = f"nodes-{_sha256(nodes)[:8]}.json"
    if not (md / nodes_file).exists():
        (md / nodes_file).write_text(json.dumps(nodes, indent=1) + "\n")

    line = {
        "run_id": manifest["run_id"],
        "timestamp": manifest["timestamp"],
        "ob_version": manifest["ob_version"],
        "benchmark": bench,
        "plan": plan_copy.name,
        "selection": selection,
        "capabilities": sorted(capabilities),
        "nodes": nodes_file,
        "node_count": len(nodes),
        "trace": trace,
    }
    with open(md / RUNLOG, "a") as fh:
        fh.write(json.dumps(line) + "\n")
    return line


def read_runs(out_dir: Path) -> List[dict]:
    """Every recorded run, oldest first.

    A line without a `benchmark` block was written by an unreleased ob in
    another format; it is skipped, not read as a run.
    """
    path = Path(out_dir) / METADATA_DIR / RUNLOG
    if not path.is_file():
        return []
    lines = (json.loads(x) for x in path.read_text().splitlines() if x.strip())
    return [run for run in lines if "benchmark" in run]


def recorded_benchmark(out_dir: Path) -> Optional[dict]:
    """The identity of the last run into *out_dir*, or None if nothing ran there.

    An out/ written before runs.jsonl existed has only the plan copy: `id` and
    `version` are read from it as plain fields; its hash is recomputed, the only
    option left, and None if the copy no longer parses.
    """
    runs = read_runs(out_dir)
    if runs:
        return runs[-1]["benchmark"]
    copy = Path(out_dir) / METADATA_DIR / PLAN_COPY
    if not copy.is_file():
        return None
    raw = yaml.safe_load(copy.read_text()) or {}
    try:
        from omnibenchmark.model import Benchmark

        summary_hash = Benchmark.from_yaml(copy).summary_hash()
    except Exception:
        summary_hash = None
    return {
        "id": raw.get("id"),
        "version": str(raw.get("version")),
        "summary_hash": summary_hash,
    }


def identity_mismatch(model, out_dir: Path, versioned: bool = False) -> Optional[str]:
    """Why *model* must not run into *out_dir*: another benchmark id or version.

    *versioned*: the storage keeps earlier versions (remote storage), so a new
    version may run into the same out dir; only the id is compared.
    """
    rec = recorded_benchmark(out_dir)
    if rec is None:
        return None
    mine = benchmark_record(model)
    diffs = [
        f"{k} {rec[k]!r} (plan: {mine[k]!r})"
        for k in (("id",) if versioned else ("id", "version"))
        if rec.get(k) != mine[k]
    ]
    if not diffs:
        return None
    return (
        f"'{out_dir}' holds results of another benchmark: {', '.join(diffs)}. "
        "Use a new --out-dir."
    )


def check_out_dir(model, out_dir: Path) -> Tuple[List[str], List[str]]:
    """(errors, notes) for reading *out_dir* with *model*.

    Errors: another benchmark id or version, or a plan edited since the last
    run. Notes: earlier runs under other plan hashes, and a partial out/.
    Nothing recorded means nothing ran: no errors, the reader reports it.
    """
    mismatch = identity_mismatch(model, out_dir)
    if mismatch:
        return [mismatch], []
    rec = recorded_benchmark(out_dir)
    if rec is None:
        return [], []
    mine = model.summary_hash()
    if rec.get("summary_hash") != mine:
        return [
            f"The plan changed since the last run into '{out_dir}' "
            f"(out/ {str(rec.get('summary_hash'))[:8]}, plan {mine[:8]}). "
            "Run `ob run` first."
        ], []

    notes = []
    runs = read_runs(out_dir)
    other = sorted({r["benchmark"]["summary_hash"][:8] for r in runs} - {mine[:8]})
    if other:
        # TODO(007 §4.6): fail on outputs orphaned by those plans, not just warn.
        notes.append(
            f"'{out_dir}' also holds runs of earlier plan revisions ({', '.join(other)}); "
            "their outputs may no longer belong to this plan."
        )
    current = [r for r in runs if r["benchmark"]["summary_hash"] == mine]
    if current and all(r["selection"] for r in current):
        picked = "; ".join(_describe(r["selection"]) for r in current)
        notes.append(f"'{out_dir}' holds a partial run of this plan: {picked}.")
    return [], notes


def refuse_foreign_out_dir(model, out_dir: Path) -> bool:
    """Log `check_out_dir`'s errors and notes; False if there were errors."""
    errors, notes = check_out_dir(model, out_dir)
    for msg in errors:
        logger.error(msg)
    if errors:
        return False
    for msg in notes:
        logger.warning(msg)
    return True


def _describe(selection: dict) -> str:
    """The CLI flag that made *selection* (see `selection_of`)."""
    if "module" in selection:
        return f"-m {selection['module']}"
    if "until" in selection:
        return f"--until {selection['until']}"
    return f"--filter {selection['filter']['sha256'][:8]}"
