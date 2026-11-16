"""Slices: ancestor-closed parts of a run, packed as tar.gz or squashfs. SPIKE.

A node key hashes everything upstream of a node (module commit, entrypoint,
environment, parameters, paths, parent keys). A slice is a set of nodes closed
under ancestors; its hash is over its node keys, so a matching slice file is
valid for the plan that computed the same hash, node for node.

`ob run` writes the keys to out/.metadata/nodes.json; packing and
`ob run --use-snapshots` read only that file.

Snapshots live in a local directory: `ob config [--local] snapshots.path DIR`,
else `snapshots/` next to the plan.
"""

import hashlib
import io
import json
import os
import shutil
import subprocess
import tarfile
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

from omnibenchmark.config import LOCAL_CONFIG
from omnibenchmark.constants import KEY_FILE

NODES_FILE = Path(".metadata") / "nodes.json"
MOUNTS_DIR = ".slices"
_APPTAINER_BIN = "/usr/libexec/apptainer/bin"  # apptainer ships squashfuse_ll


def store_dir(plan_dir: Path) -> Path:
    """The snapshots directory: the plan's .ob/config, then the user's
    config, then `snapshots/`. A relative path is relative to the plan."""
    from omnibenchmark.config import ConfigAccessor

    for cfg in (ConfigAccessor(plan_dir / LOCAL_CONFIG), ConfigAccessor()):
        path = cfg.get("snapshots", "path")
        if path:
            return plan_dir / Path(path).expanduser()
    return plan_dir / "snapshots"


def _sha(obj) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True).encode()).hexdigest()


def _parent_ids(node) -> List[str]:
    ids = set(node.parents) | set(node.gathered_from)
    if node.parent_id:
        ids.add(node.parent_id)
    return sorted(ids)


def _env_digest(module, out_dir: Path) -> Optional[str]:
    env = module.resolved_environment
    if env is None or not env.reference:
        return None
    path = out_dir / env.reference
    if path.is_file():  # conda yaml: hash the contents, the path is local
        return hashlib.sha256(path.read_bytes()).hexdigest()
    return f"{env.backend_type}:{env.reference}"


def node_keys(nodes, out_dir: Path) -> Dict[str, str]:
    """{node id: key}. A key commits to the node's whole upstream cone."""
    by_id = {n.id: n for n in nodes}
    keys: Dict[str, str] = {}

    def key(n) -> str:
        if n.id not in keys:
            m = n.module
            keys[n.id] = _sha(
                {
                    "repo": m.repository_url,
                    # a working copy has no commit to pin: never reusable
                    "commit": m.commit if not m.dirty else os.urandom(8).hex(),
                    "entrypoint": str(m.entrypoint),
                    "env": _env_digest(m, out_dir),
                    "params": n.parameters.serialize() if n.parameters else None,
                    "dir": n.node_dir,
                    "outputs": n.outputs,
                    "parents": sorted(
                        key(by_id[p]) for p in _parent_ids(n) if p in by_id
                    ),
                }
            )[:16]
        return keys[n.id]

    for n in nodes:
        key(n)
    return keys


def write_nodes(nodes, model, out_dir: Path) -> None:
    keys = node_keys(nodes, out_dir)
    meta = {
        "benchmark": model.id,
        "version": model.get_version(),
        "plan_hash": model.summary_hash(),
        "stages": list(model.get_stages()),
        "nodes": [
            {
                "id": n.id,
                "stage": n.stage_id,
                "dir": n.node_dir,
                "key": keys[n.id],
                "outputs": sorted(n.outputs.values()),
                "parents": _parent_ids(n),
            }
            for n in nodes
        ],
    }
    path = out_dir / NODES_FILE
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(meta, indent=1))


def read_nodes(out_dir: Path) -> dict:
    return json.loads((out_dir / NODES_FILE).read_text())


def _ancestors(meta: dict) -> Dict[str, set]:
    """{node id: ids of all its ancestors}."""
    by_id = {n["id"]: n for n in meta["nodes"]}
    done: Dict[str, set] = {}

    def of(node_id: str) -> set:
        if node_id not in done:
            done[node_id] = set()  # a plan is acyclic; this only guards a bad file
            parents = [p for p in by_id[node_id]["parents"] if p in by_id]
            done[node_id] = set(parents).union(*(of(p) for p in parents))
        return done[node_id]

    for n in meta["nodes"]:
        of(n["id"])
    return done


def select(
    meta: dict,
    until: Optional[str] = None,
    root: Optional[str] = None,
    ancestors: Optional[Dict[str, set]] = None,
) -> List[dict]:
    """Nodes up to stage `until`, in the group of the node `root`, plus their
    ancestors.

    A node is in `root`'s group when `root` is its only ancestor at `root`'s
    stage. A gather over several groups is therefore in none of them.
    """
    order = {s: i for i, s in enumerate(meta["stages"])}
    by_id = {n["id"]: n for n in meta["nodes"]}
    anc = ancestors if ancestors is not None else _ancestors(meta)

    def in_group(n: dict) -> bool:
        if root is None or n["id"] == root:
            return True
        stage = by_id[root]["stage"]
        return {a for a in anc[n["id"]] if by_id[a]["stage"] == stage} == {root}

    picked = {
        n["id"]
        for n in meta["nodes"]
        if (until is None or order.get(n["stage"], len(order)) <= order[until])
        and in_group(n)
    }
    for node_id in list(picked):
        picked |= anc[node_id]
    return [n for n in meta["nodes"] if n["id"] in picked]


def slice_hash(nodes: List[dict]) -> str:
    return _sha(sorted(n["key"] for n in nodes))[:8]


def slice_label(until: Optional[str], root: Optional[dict]) -> str:
    """`DATA.d1.94686c86.until-PCA`: the root's directory segment, then the cut."""
    segment = os.path.basename(root["dir"]) if root else None
    return ".".join(filter(None, [segment, until and f"until-{until}"])) or "all"


def _node_files(out_dir: Path, nodes: List[dict], all_dirs: set) -> List[str]:
    """Files of the given nodes, relative to out_dir. A child node's directory
    sits inside its parent's, so other nodes' dirs are pruned from the walk."""
    files = []
    for n in nodes:
        for root, dirs, names in os.walk(out_dir / n["dir"]):
            rel = os.path.relpath(root, out_dir)
            dirs[:] = [d for d in dirs if os.path.join(rel, d) not in all_dirs]
            files += [os.path.join(rel, name) for name in names]
    return sorted(files)


def _check_produced(out_dir: Path, nodes: List[dict]) -> None:
    """Every node's directory must hold the outputs the plan declares and the
    key the plan computes now. The key is left there by the rule that ran, so
    a directory produced under another commit, environment, parameter or
    parent is caught, as is one that never ran."""
    stale, incomplete = [], []
    for n in nodes:
        key_file = out_dir / n["dir"] / KEY_FILE
        if not key_file.is_file() or key_file.read_text().strip() != n["key"]:
            stale.append(n["dir"])
        elif not all((out_dir / o).exists() for o in n.get("outputs", [])):
            incomplete.append(n["dir"])
    for dirs, why in (
        (stale, "not produced by the plan as it resolves now"),
        (incomplete, "missing declared outputs"),
    ):
        if dirs:
            raise ValueError(
                f"{len(dirs)} node(s) {why}, e.g. {sorted(dirs)[0]}. "
                "Run the plan to completion, then archive."
            )


def _relative(value, prefix: str):
    """`value` with every string under `prefix` made relative to it."""
    if isinstance(value, str):
        return value.removeprefix(prefix)
    if isinstance(value, list):
        return [_relative(v, prefix) for v in value]
    if isinstance(value, dict):
        return {k: _relative(v, prefix) for k, v in value.items()}
    return value


def _add(tar: tarfile.TarFile, out_dir: Path, f: str) -> None:
    """Add out_dir/f as f. A lineage.json names its member files by absolute
    path, for the module that reads it; in a snapshot they are relative to
    out/, so the record holds on another machine."""
    if os.path.basename(f) != "lineage.json":
        tar.add(out_dir / f, arcname=f)
        return
    prefix = str(out_dir.resolve()) + os.sep
    body = json.loads((out_dir / f).read_text())
    data = json.dumps(_relative(body, prefix), indent=1).encode()
    info = tar.gettarinfo(out_dir / f, arcname=f)
    info.size = len(data)
    tar.addfile(info, io.BytesIO(data))


def pack(
    out_dir: Path,
    dest: Path,
    until: Optional[str] = None,
    group_by: Optional[str] = None,
    squashfs: bool = False,
    base_until: Optional[str] = None,
) -> List[str]:
    """Pack one snapshot, or with `group_by` one per node of that stage.
    Returns one line per snapshot. A snapshot that is already there is left
    alone: the same name is the same content."""
    meta = read_nodes(out_dir)
    anc = _ancestors(meta)
    roots: List[Optional[dict]] = [None]
    if group_by:
        roots = [n for n in meta["nodes"] if n["stage"] == group_by]
        if not roots:
            raise ValueError(f"no node of stage {group_by!r} to group by")
    report = []
    for root in roots:
        target = _pack_one(meta, anc, out_dir, dest, until, root, squashfs, base_until)
        report.append(target)
    return report


def _pack_one(meta, anc, out_dir, dest, until, root, squashfs, base_until) -> str:
    """With `base_until`, the file holds only the nodes above that cut and is
    named `..._<hash>~<base hash>`: the hash is still the whole slice's, and
    the base snapshot has to be used with it."""
    root_id = root["id"] if root else None
    nodes = select(meta, until, root_id, anc)
    _check_produced(out_dir, nodes)
    name = (
        f"{meta['benchmark']}-{meta['version']}-{meta['plan_hash'][:8]}"
        f"-{slice_label(until, root)}_{slice_hash(nodes)}"
    )
    whole, needs = len(nodes), ""
    if base_until:
        # below the group's own stage the base is the ungrouped slice
        base = select(meta, base_until, root_id, anc) or []
        if not any(n["id"] == root_id for n in base):
            base = select(meta, base_until, None, anc)
        needs = slice_hash(base)
        name += f"~{needs}"
        nodes = [n for n in nodes if n not in base]
    if not nodes:
        raise ValueError("the slice selects no nodes")
    files = _node_files(out_dir, nodes, {n["dir"] for n in meta["nodes"]})
    dest.mkdir(parents=True, exist_ok=True)
    target = dest / (name + (".sqfs" if squashfs else ".tar.gz"))
    if target.exists():
        return f"exists  {target}"
    part = target.with_name(target.name + ".part")
    part.unlink(missing_ok=True)
    try:
        if squashfs:
            # mksquashfs >= 4.5 reads a tar stream, so both formats share the
            # member list below
            proc = subprocess.Popen(
                ["mksquashfs", "-", str(part), "-tar", "-comp", "zstd", "-quiet"]
                + ["-no-progress"],
                stdin=subprocess.PIPE,
                stdout=subprocess.DEVNULL,
            )
            with tarfile.open(fileobj=proc.stdin, mode="w|") as tar:
                for f in files:
                    _add(tar, out_dir, f)
            assert proc.stdin is not None
            proc.stdin.close()
            if proc.wait() != 0:
                raise RuntimeError("mksquashfs failed")
        else:
            with tarfile.open(part, "w:gz") as tar:
                for f in files:
                    _add(tar, out_dir, f)
        os.replace(part, target)
    finally:
        part.unlink(missing_ok=True)
    _index(meta, target, nodes, whole, root, until, needs)
    return f"created {target}"


INDEX = "README.md"
_INDEX_HEAD = """# Snapshots

One row per file. A file matches a plan when the plan computes the same slice
hash (the 8 characters after the last `_`). A file that names a base holds
only what the base lacks, and is used together with it.

| File | Benchmark | Slice | Holds | Size | Needs | Packed |
|---|---|---|---|---|---|---|
"""


def _index(meta, target: Path, nodes, whole: int, root, until, needs: str) -> None:
    """Append the new snapshot to the store's README.md, for people. Nothing
    reads it back."""
    stages = [s for s in meta["stages"] if any(n["stage"] == s for n in nodes)]
    held = f"{len(nodes)} of the slice's {whole} nodes" if needs else f"{whole} nodes"
    what = " and ".join(
        filter(
            None,
            [
                f"below `{os.path.basename(root['dir'])}`" if root else "",
                f"up to `{until}`" if until else "",
            ],
        )
    )
    row = [
        f"`{target.name}`",
        f"{meta['benchmark']} {meta['version']} (plan {meta['plan_hash'][:8]})",
        what or "the whole run",
        f"{held}; stages {', '.join(stages)}",
        f"{target.stat().st_size / 1e6:.1f} MB",
        f"the file ending `_{needs}`" if needs else "",
        datetime.now().astimezone().isoformat(timespec="seconds"),
    ]
    index = target.parent / INDEX
    if not index.exists():
        index.write_text(_INDEX_HEAD)
    with index.open("a") as f:
        f.write("| " + " | ".join(row) + " |\n")


def _squashfuse() -> Optional[str]:
    for name in ("squashfuse_ll", "squashfuse"):
        found = shutil.which(name) or shutil.which(name, path=_APPTAINER_BIN)
        if found and os.path.exists("/dev/fuse"):
            return found
    return None


def _mount(image: Path, out_dir: Path) -> bool:
    """Mount `image` read-only under out/.slices and symlink its files into
    out/. Directories are real, so new child nodes can be written next to the
    cached files; the cached files themselves cannot be written."""
    tool = _squashfuse()
    if tool is None:
        return False
    mnt = out_dir / MOUNTS_DIR / image.name
    mnt.mkdir(parents=True, exist_ok=True)
    if not os.path.ismount(mnt):
        if subprocess.run([tool, str(image.resolve()), str(mnt)]).returncode != 0:
            return False
    for root, _dirs, names in os.walk(mnt):
        rel = os.path.relpath(root, mnt)
        (out_dir / rel).mkdir(parents=True, exist_ok=True)
        for name in names:
            link = out_dir / rel / name
            if not link.exists() and not link.is_symlink():
                src = os.path.join(root, name)
                link.symlink_to(os.path.relpath(src, link.parent))
                # the link carries the file's mtime, or every cached file would
                # look newer than whatever is computed next
                st = os.stat(src)
                times = (st.st_atime_ns, st.st_mtime_ns)
                os.utime(link, ns=times, follow_symlinks=False)
    return True


def _extract(image: Path, out_dir: Path) -> None:
    if image.name.endswith(".sqfs"):
        # no fuse: unsquashfs needs none
        subprocess.run(
            ["unsquashfs", "-f", "-quiet", "-no-progress", "-d", str(out_dir)]
            + [str(image)],
            check=True,
            stdout=subprocess.DEVNULL,
        )
    else:
        with tarfile.open(image) as tar:
            tar.extractall(out_dir, filter="data")


def stored_keys(image: Path, dirs: List[str]) -> Dict[str, str]:
    """{node dir: the key its rule left in the snapshot}, for the dirs that
    have one. Read from the file itself: no mount, nothing extracted."""
    keys: Dict[str, str] = {}
    if image.name.endswith(".sqfs"):
        for d in dirs:  # one call each: a missing file fails only its own
            cat = subprocess.run(
                ["unsquashfs", "-cat", str(image), f"{d}/{KEY_FILE}"],
                capture_output=True,
                text=True,
            )
            if cat.returncode == 0:
                keys[d] = cat.stdout.strip()
        return keys
    wanted = {f"{d}/{KEY_FILE}": d for d in dirs}
    with tarfile.open(image) as tar:
        for member in tar:
            if member.name in wanted:
                data = tar.extractfile(member)
                keys[wanted[member.name]] = data.read().decode().strip() if data else ""
    return keys


def use_snapshots(out_dir: Path, source: Path, apply: bool = True) -> List[str]:
    """Fill out_dir from every snapshot in `source` whose slice hash this
    plan computes too: a squashfs image is mounted read-only when it can be,
    anything else is extracted. Returns one line per snapshot; with
    `apply=False` it only says what would be used."""
    meta = read_nodes(out_dir)
    anc = _ancestors(meta)
    # A grouped snapshot carries its root's directory segment in its name, so
    # only the nodes some file names can be roots.
    names = [p.name for p in source.glob("*") if p.is_file()]
    roots = [
        n["id"]
        for n in meta["nodes"]
        if any(f"-{os.path.basename(n['dir'])}" in name for name in names)
    ]
    wanted: Dict[str, List[dict]] = {}
    for root in [None, *roots]:
        for until in [None, *meta["stages"]]:
            nodes = select(meta, until, root, anc)
            if nodes:
                wanted.setdefault(slice_hash(nodes), nodes)

    # hash -> (file, base hash or ""); a file named `_<hash>~<base>` holds
    # only what `<base>` lacks
    found: Dict[str, tuple] = {}
    exts = (".sqfs", ".tar.gz") if _squashfuse() else (".tar.gz", ".sqfs")
    for ext in reversed(exts):  # the preferred format overwrites the other
        for path in sorted(source.glob(f"*{ext}"), key=lambda p: "~" in p.name):
            h, _, base = path.name[: -len(ext)].rsplit("_", 1)[-1].partition("~")
            if h in wanted and (not base or base in wanted):
                found[h] = (path, base)

    report, covered = [], set()
    for h in sorted(found, key=lambda h: -len(wanted[h])):  # biggest first
        image, base = found[h]
        ids = {n["id"] for n in wanted[h]}
        if base:
            if base not in found:
                report.append(f"skipped {image.name}: its base slice is missing")
                continue
            ids -= {n["id"] for n in wanted[base]}
        if ids <= covered:
            continue
        how = "extracted"
        if not apply:
            how = "would use"
            covered |= ids
            report.append(f"{how} {image.name}: {len(ids)} nodes")
            continue
        # The name only promises a slice hash. Each node's own key, left in
        # its directory by the rule that ran, has to be the one this plan
        # computes: that is what makes a node reusable, not the plan's hash.
        expected = {n["dir"]: n["key"] for n in wanted[h] if n["id"] in ids}
        stored = stored_keys(image, sorted(expected))
        wrong = sorted(d for d, key in expected.items() if stored.get(d) != key)
        if wrong:
            report.append(
                f"skipped {image.name}: {len(wrong)} of its {len(expected)} nodes "
                f"do not carry this plan's key, e.g. {wrong[0]}"
            )
            continue
        if image.name.endswith(".sqfs") and _mount(image, out_dir):
            how = "mounted"
        else:
            _extract(image, out_dir)
        covered |= ids
        report.append(f"{how} {image.name}: {len(ids)} nodes")
    return report


def unmount(out_dir: Path) -> List[str]:
    done = []
    for mnt in sorted((out_dir / MOUNTS_DIR).glob("*")):
        if os.path.ismount(mnt):
            # busy (a shell sitting in it, another mount namespace): detach
            # now, the kernel finishes when the last user leaves
            if subprocess.run(["fusermount", "-u", str(mnt)]).returncode != 0:
                subprocess.run(["fusermount", "-uz", str(mnt)], check=True)
            done.append(mnt.name)
    return done
