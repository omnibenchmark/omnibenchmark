"""Slices: ancestor-closed parts of a run, packed as tar.gz or squashfs. SPIKE.

A node key hashes everything upstream of a node (module commit, entrypoint,
environment, parameters, paths, parent keys). A slice is a set of nodes closed
under ancestors; its hash is over its node keys, so a matching slice file is
valid for the plan that computed the same hash, node for node.

`ob run` writes the keys to out/.metadata/nodes.json. `pack` and `seed` read
only that file, so neither resolves a plan.
"""

import hashlib
import json
import os
import shutil
import subprocess
import tarfile
from pathlib import Path
from typing import Dict, List, Optional

NODES_FILE = Path(".metadata") / "nodes.json"
MOUNTS_DIR = ".slices"
_APPTAINER_BIN = "/usr/libexec/apptainer/bin"  # apptainer ships squashfuse_ll


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
                "dataset": (
                    n.template_context.provides.get("dataset")
                    if n.template_context
                    else None
                ),
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


def select(meta: dict, until: Optional[str], dataset: Optional[str]) -> List[dict]:
    """Nodes up to stage `until` and/or in `dataset`'s subtree, plus ancestors."""
    order = {s: i for i, s in enumerate(meta["stages"])}
    by_id = {n["id"]: n for n in meta["nodes"]}
    todo = [
        n
        for n in meta["nodes"]
        if (until is None or order.get(n["stage"], len(order)) <= order[until])
        and (dataset is None or n["dataset"] == dataset)
    ]
    picked: Dict[str, dict] = {}
    while todo:
        n = todo.pop()
        if n["id"] not in picked:
            picked[n["id"]] = n
            todo += [by_id[p] for p in n["parents"] if p in by_id]
    return list(picked.values())


def slice_hash(nodes: List[dict]) -> str:
    return _sha(sorted(n["key"] for n in nodes))[:8]


def slice_label(until: Optional[str], dataset: Optional[str]) -> str:
    return ".".join(filter(None, [dataset, until and f"until-{until}"])) or "all"


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


def pack(
    out_dir: Path,
    dest: Path,
    until: Optional[str] = None,
    dataset: Optional[str] = None,
    squashfs: bool = False,
) -> Path:
    meta = read_nodes(out_dir)
    nodes = select(meta, until, dataset)
    if not nodes:
        raise ValueError("the slice selects no nodes")
    files = _node_files(out_dir, nodes, {n["dir"] for n in meta["nodes"]})
    name = (
        f"{meta['benchmark']}-{meta['version']}-{meta['plan_hash'][:8]}"
        f"-{slice_label(until, dataset)}_{slice_hash(nodes)}"
    )
    dest.mkdir(parents=True, exist_ok=True)
    target = dest / (name + (".sqfs" if squashfs else ".tar.gz"))
    if target.exists():
        raise FileExistsError(f"{target} exists; a slice is written once")
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
                    tar.add(out_dir / f, arcname=f)
            assert proc.stdin is not None
            proc.stdin.close()
            if proc.wait() != 0:
                raise RuntimeError("mksquashfs failed")
        else:
            with tarfile.open(part, "w:gz") as tar:
                for f in files:
                    tar.add(out_dir / f, arcname=f)
        os.replace(part, target)
    finally:
        part.unlink(missing_ok=True)
    return target


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


def seed(out_dir: Path, source: Path, squashfs: bool = False) -> List[str]:
    """Seed out_dir from every slice file in `source` whose hash this plan
    computes too. Returns one line per slice used."""
    meta = read_nodes(out_dir)
    datasets = sorted({n["dataset"] for n in meta["nodes"] if n["dataset"]})
    wanted: Dict[str, List[dict]] = {}
    for dataset in [None, *datasets]:
        for until in [None, *meta["stages"]]:
            nodes = select(meta, until, dataset)
            if nodes:
                wanted.setdefault(slice_hash(nodes), nodes)

    found: Dict[str, Path] = {}
    exts = (".sqfs", ".tar.gz") if squashfs else (".tar.gz",)
    for ext in reversed(exts):  # the preferred format overwrites the other
        for path in source.glob(f"*_????????{ext}"):
            h = path.name[: -len(ext)].rsplit("_", 1)[1]
            if h in wanted:
                found[h] = path

    report, covered = [], set()
    for h in sorted(found, key=lambda h: -len(wanted[h])):  # biggest first
        ids = {n["id"] for n in wanted[h]}
        if ids <= covered:
            continue
        image = found[h]
        how = "extracted"
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
