import pytest

from omnibenchmark.slices import select, slice_hash, slice_label


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
