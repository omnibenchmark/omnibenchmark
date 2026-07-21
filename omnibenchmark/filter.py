"""obfilter — apply obeditor-style picks to prune the benchmark DAG.

SPIKE (see scratch/obfilter_native_plan.md). Picks shape (blob v3):

    {stage_id: {(module_id | "*"): "all" | "first" | ["<combo_hash8>", ...]}}

- A stage absent from picks is dropped entirely.
- Within a selected stage, an explicit module key wins; otherwise "*" (if present)
  applies to any module. "*" cannot be mixed with explicit module keys (raises).
- spec: "all" keep every combo · "first" keep first combo · [hashes] keep those combos
  (by `Params.hash_short()`).

The blob anchor `parent.sha256` is `Benchmark.summary_hash()` (execution-only), so
administrative edits to the parent don't break a slice.
"""

import base64
import gzip
import json

BLOB_VERSION = 3
WILDCARD = "*"
DROP = None  # module_spec() returns this when a (stage, module) is not selected


class FilterError(ValueError):
    """Malformed picks, or an unsupported/undecodable blob."""


# --------------------------------------------------------------------------- codec


def pack_blob(picks, parent):
    """picks + parent -> packed transport string (canonical JSON → gzip → urlsafe b64, no pad)."""
    validate_picks(picks)
    blob = {"v": BLOB_VERSION, "parent": parent, "picks": picks}
    raw = json.dumps(blob, sort_keys=True, separators=(",", ":")).encode("utf-8")
    packed = base64.urlsafe_b64encode(gzip.compress(raw, compresslevel=9, mtime=0))
    return packed.decode("ascii").rstrip("=")


def unpack_blob(packed):
    """Inverse of pack_blob. Returns the blob dict; raises FilterError on bad input/version."""
    packed = packed.strip()
    pad = "=" * (-len(packed) % 4)
    try:
        raw = gzip.decompress(base64.urlsafe_b64decode(packed + pad))
        blob = json.loads(raw)
    except Exception as e:
        raise FilterError(f"could not decode filter blob: {e}") from e
    if blob.get("v") != BLOB_VERSION:
        raise FilterError(
            f"unsupported filter blob version {blob.get('v')!r} (expected {BLOB_VERSION})"
        )
    validate_picks(blob.get("picks") or {})
    return blob


# --------------------------------------------------------------------------- schema


def validate_picks(picks):
    """Structural validation. Raises FilterError; '*' cannot mix with explicit keys."""
    if not isinstance(picks, dict):
        raise FilterError("picks must be a mapping stage_id -> {module_id|'*': spec}")
    for stage_id, mods in picks.items():
        if not isinstance(mods, dict):
            raise FilterError(
                f"picks[{stage_id!r}] must be a mapping module_id|'*' -> spec"
            )
        if WILDCARD in mods and len(mods) > 1:
            raise FilterError(
                f"stage {stage_id!r}: '*' wildcard cannot be mixed with explicit module keys"
            )
        for key, spec in mods.items():
            if spec in ("all", "first"):
                continue
            if isinstance(spec, list) and all(isinstance(h, str) for h in spec):
                continue
            raise FilterError(
                f"picks[{stage_id!r}][{key!r}]: spec must be 'all', 'first', or a list of hashes"
            )


# --------------------------------------------------------------------------- apply


def module_spec(picks, stage_id, module_id):
    """Pick spec for (stage, module), or DROP (None) if not selected.

    Explicit module key wins; else the stage's "*" wildcard; else DROP.
    """
    stage = picks.get(stage_id)
    if not stage:
        return DROP
    if module_id in stage:
        return stage[module_id]
    if WILDCARD in stage:
        return stage[WILDCARD]
    return DROP


def keeps_module(picks, stage_id, module_id):
    return module_spec(picks, stage_id, module_id) is not DROP


def filter_params(params_list, spec):
    """Apply a combo spec to an expanded params list.

    "all"/no-params -> unchanged · "first" -> first combo · [hashes] -> matching combos.
    """
    if spec == "all" or params_list == [None]:
        return params_list
    if spec == "first":
        return params_list[:1]
    wanted = set(spec)
    return [p for p in params_list if p.hash_short() in wanted]


# --------------------------------------------------------------------------- drift


def find_orphans(picks, model):
    """Human-readable descriptions of picks that don't resolve against `model`.

    Renamed/removed stage or module, or an explicit combo hash that no longer expands.
    "*" and "all"/"first" specs never orphan.
    """
    from omnibenchmark.model.params import Params

    orphans = []
    stages_by_id = {s.id: s for s in model.stages}
    for stage_id, mods in picks.items():
        stage = stages_by_id.get(stage_id)
        if stage is None:
            orphans.append(f"stage '{stage_id}' (not in benchmark)")
            continue
        modules_by_id = {m.id: m for m in stage.modules}
        for key, spec in mods.items():
            if key == WILDCARD:
                continue
            module = modules_by_id.get(key)
            if module is None:
                orphans.append(f"module '{stage_id}/{key}' (not in stage)")
                continue
            if isinstance(spec, list):
                have = set()
                for pset in module.parameters or []:
                    for p in Params.expand_from_parameter(pset):
                        have.add(p.hash_short())
                for h in spec:
                    if h not in have:
                        orphans.append(
                            f"combo '{stage_id}/{key}/{h}' (no longer expands)"
                        )
    return orphans
