import io
import json
import tarfile

import pytest

from omnibenchmark.slices import (
    NODES_FILE,
    select,
    slice_hash,
    slice_label,
    use_snapshots,
)


def _meta():
    def n(id, stage, parents, key=None):
        return dict(
            id=id, stage=stage, dir=f"x/{stage}.{id}", key=key or id, parents=parents
        )

    return {
        "stages": ["fetch", "data", "method", "summary"],
        "nodes": [
            n("f", "fetch", []),
            n("d1", "data", ["f"]),
            n("d2", "data", ["f"]),
            n("m1", "method", ["d1"]),
            n("m2", "method", ["d2"]),
            n("s", "summary", ["m1", "m2"]),
        ],
    }


def _ids(nodes):
    return sorted(n["id"] for n in nodes)


@pytest.mark.short
def test_select_is_ancestor_closed():
    meta = _meta()
    assert _ids(select(meta, "data")) == ["d1", "d2", "f"]
    # the trunk comes along; the gather over both groups is in neither
    assert _ids(select(meta, None, "d1")) == ["d1", "f", "m1"]
    assert _ids(select(meta, "data", "d1")) == ["d1", "f"]
    # a cut below the root's stage holds nothing of its group
    assert select(meta, "fetch", "d1") == []
    assert len(select(meta)) == 6
    d1 = meta["nodes"][1]
    assert slice_label("method", d1) == "data.d1.until-method"
    assert slice_label(None, None) == "all"


@pytest.mark.short
def test_slice_hash_follows_keys_not_order():
    meta = _meta()
    cut = slice_hash(select(meta, "data"))
    meta["nodes"].reverse()
    assert slice_hash(select(meta, "data")) == cut
    meta["nodes"][0]["key"] = "changed"  # the summary node: outside the cut
    assert slice_hash(select(meta, "data")) == cut
    next(n for n in meta["nodes"] if n["id"] == "f")["key"] = "changed"
    assert slice_hash(select(meta, "data")) != cut


def _tar(path, members: dict):
    with tarfile.open(path, "w:gz") as tar:
        for name, text in members.items():
            info = tarfile.TarInfo(name)
            info.size = len(text)
            tar.addfile(info, io.BytesIO(text.encode()))


@pytest.mark.short
def test_use_snapshots_checks_each_node_key(tmp_path):
    """A file named by the right slice hash is used only if every node in it
    carries the key this plan computes."""
    meta = {
        "stages": ["data"],
        "nodes": [
            dict(id="d1", stage="data", dir="data.d1", key="k1", parents=[]),
            dict(id="d2", stage="data", dir="data.d2", key="k2", parents=[]),
        ],
    }
    out, store = tmp_path / "out", tmp_path / "store"
    (out / NODES_FILE).parent.mkdir(parents=True)
    (out / NODES_FILE).write_text(json.dumps(meta))
    store.mkdir()
    name = f"b-1-00000000-all_{slice_hash(meta['nodes'])}.tar.gz"
    good = {"data.d1/.key": "k1\n", "data.d1/x": "1", "data.d2/.key": "k2\n"}

    _tar(store / name, {**good, "data.d2/.key": "other\n"})
    (line,) = use_snapshots(out, store)
    assert line.startswith("skipped") and "data.d2" in line
    assert not (out / "data.d1").exists()

    _tar(store / name, {k: v for k, v in good.items() if k != "data.d2/.key"})
    assert use_snapshots(out, store)[0].startswith("skipped")

    _tar(store / name, good)
    assert use_snapshots(out, store) == [f"extracted {name}: 2 nodes"]
    assert (out / "data.d1" / "x").read_text() == "1"
