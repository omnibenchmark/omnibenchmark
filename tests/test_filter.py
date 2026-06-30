"""Unit tests for the obfilter spike (omnibenchmark/filter.py)."""

import pytest

from omnibenchmark import filter as f
from omnibenchmark.model.params import Params


# --------------------------------------------------------------------------- codec


@pytest.mark.short
class TestBlobCodec:
    def test_round_trip(self):
        picks = {"data": {"*": "all"}, "methods": {"M1": "all", "M2": ["abc12345"]}}
        parent = {"sha256": "de" * 32, "url": "https://x/b"}
        blob = f.unpack_blob(f.pack_blob(picks, parent))
        assert blob["v"] == 3
        assert blob["picks"] == picks
        assert blob["parent"] == parent

    def test_pack_is_deterministic(self):
        picks = {"s": {"m": "all"}}
        parent = {"sha256": "00" * 32}
        assert f.pack_blob(picks, parent) == f.pack_blob(picks, parent)

    def test_unpack_rejects_garbage(self):
        with pytest.raises(f.FilterError):
            f.unpack_blob("not-a-real-blob!!")

    def test_unpack_rejects_wrong_version(self):
        import base64
        import gzip
        import json

        raw = json.dumps({"v": 2, "parent": {}, "picks": {}}).encode()
        packed = (
            base64.urlsafe_b64encode(gzip.compress(raw, mtime=0)).decode().rstrip("=")
        )
        with pytest.raises(f.FilterError, match="version"):
            f.unpack_blob(packed)


# --------------------------------------------------------------------------- schema


@pytest.mark.short
class TestValidatePicks:
    def test_accepts_specs(self):
        f.validate_picks({"s": {"m1": "all", "m2": "first", "m3": ["a", "b"]}})

    def test_wildcard_alone_ok(self):
        f.validate_picks({"s": {"*": "all"}})

    def test_wildcard_mixed_raises(self):
        with pytest.raises(f.FilterError, match="wildcard cannot be mixed"):
            f.validate_picks({"s": {"*": "all", "m1": "first"}})

    def test_bad_spec_raises(self):
        with pytest.raises(f.FilterError, match="spec must be"):
            f.validate_picks({"s": {"m": "second"}})

    def test_non_mapping_raises(self):
        with pytest.raises(f.FilterError):
            f.validate_picks({"s": ["m1", "m2"]})


# --------------------------------------------------------------------------- apply


@pytest.mark.short
class TestModuleSpec:
    PICKS = {"data": {"*": "first"}, "methods": {"M1": "all"}}

    def test_explicit_module(self):
        assert f.module_spec(self.PICKS, "methods", "M1") == "all"

    def test_wildcard_matches_any(self):
        assert f.module_spec(self.PICKS, "data", "whatever") == "first"

    def test_unpicked_module_drops(self):
        assert f.module_spec(self.PICKS, "methods", "M2") is f.DROP

    def test_unpicked_stage_drops(self):
        assert f.module_spec(self.PICKS, "ghost", "X") is f.DROP

    def test_keeps_module(self):
        assert f.keeps_module(self.PICKS, "methods", "M1") is True
        assert f.keeps_module(self.PICKS, "methods", "M2") is False


@pytest.mark.short
class TestFilterParams:
    def _params(self):
        return [Params({"n": str(i)}) for i in range(3)]

    def test_all_unchanged(self):
        ps = self._params()
        assert f.filter_params(ps, "all") == ps

    def test_no_params_unchanged(self):
        assert f.filter_params([None], "first") == [None]

    def test_first_keeps_one(self):
        ps = self._params()
        assert f.filter_params(ps, "first") == ps[:1]

    def test_hash_list_keeps_matching(self):
        ps = self._params()
        wanted = ps[2].hash_short()
        kept = f.filter_params(ps, [wanted])
        assert [p.hash_short() for p in kept] == [wanted]

    def test_hash_list_no_match_empty(self):
        assert f.filter_params(self._params(), ["deadbeef"]) == []


# --------------------------------------------------------------------------- drift


class _StubParam:
    """Minimal Parameter stand-in for Params.expand_from_parameter."""

    def __init__(self, values):
        self.values = values


@pytest.mark.short
class TestFindOrphans:
    def _model(self):
        # Build a tiny stand-in model: stages with .id and .modules (.id, .parameters)
        from types import SimpleNamespace as NS

        m1 = NS(id="M1", parameters=None)
        m2 = NS(id="M2", parameters=None)
        data = NS(id="data", modules=[NS(id="D1", parameters=None)])
        methods = NS(id="methods", modules=[m1, m2])
        return NS(stages=[data, methods])

    def test_no_orphans_when_all_resolve(self):
        picks = {"data": {"*": "all"}, "methods": {"M1": "all"}}
        assert f.find_orphans(picks, self._model()) == []

    def test_missing_stage(self):
        orphans = f.find_orphans({"ghost": {"X": "all"}}, self._model())
        assert any("stage 'ghost'" in o for o in orphans)

    def test_missing_module(self):
        orphans = f.find_orphans({"methods": {"M9": "all"}}, self._model())
        assert any("methods/M9" in o for o in orphans)

    def test_wildcard_never_orphans(self):
        assert f.find_orphans({"methods": {"*": "all"}}, self._model()) == []
