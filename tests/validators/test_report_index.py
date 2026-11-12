"""Unit test for the validation read side (the status report reads it).

A verdict is keyed by the resolved file, so a path through a symlink (the
`out/by-params/` view, or the pre-0.8 nested layout) and the real path find the
same verdict.
"""

import json
import os

import pytest

from omnibenchmark.core.validators import load_verdicts, verdict_for


@pytest.mark.short
def test_validation_index_bridges_symlink_and_hash(tmp_path):
    out = tmp_path / "out"
    real_dir = out / "data" / "D2" / ".b43a8b63"  # hidden hash dir (real)
    real_dir.mkdir(parents=True)
    real_file = real_dir / "D2_data.json"
    real_file.write_text("x")
    # human-readable symlink dir -> hash dir (as omnibenchmark lays it out)
    symlink_dir = out / "data" / "D2" / "evaluate-2+2"
    symlink_dir.symlink_to(real_dir, target_is_directory=True)

    # a result keyed by the symlink path, as `ob validate outputs` records it
    vdir = out / ".validation" / "data" / "D2" / "evaluate-2+2"
    vdir.mkdir(parents=True)
    (vdir / "D2_data.json.json").write_text(
        json.dumps(
            {
                "output_file": "data/D2/evaluate-2+2/D2_data.json",
                "exit_code": 0,
                "log": ".validation/data/D2/evaluate-2+2/D2_data.json.log",
            }
        )
    )

    index = load_verdicts(out)

    # both aliases of the same real file resolve to the recorded verdict
    assert verdict_for(real_file, index)[0] == "pass"
    assert verdict_for(symlink_dir / "D2_data.json", index)[0] == "pass"
    # the log path is carried through (absolute)
    assert verdict_for(real_file, index)[1].endswith("D2_data.json.log")
    # an output with no result -> not_validated
    assert verdict_for(out / "data" / "D1" / "x.json", index)[0] == "not_validated"


def _verdict(out, rel, exit_code, i=0):
    result = out / ".validation" / f"{rel}.{i}.json"
    result.parent.mkdir(parents=True, exist_ok=True)
    result.write_text(json.dumps({"output_file": rel, "exit_code": exit_code}))
    return result


@pytest.mark.short
def test_validation_index_stale_missing_and_error(tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    (out / "a.txt").write_text("x")
    (out / "b.txt").write_text("x")
    result = _verdict(out, "a.txt", 0)
    _verdict(out, "b.txt", 127)
    _verdict(out, "gone.txt", 0)

    index = load_verdicts(out)
    assert verdict_for(out / "a.txt", index)[0] == "pass"
    assert verdict_for(out / "b.txt", index)[0] == "error"
    # a deleted output has no verdict
    assert verdict_for(out / "gone.txt", index)[0] == "not_validated"

    # rewriting the output after validation makes the verdict stale
    t = result.stat().st_mtime
    os.utime(out / "a.txt", (t + 10, t + 10))
    index = load_verdicts(out)
    assert verdict_for(out / "a.txt", index)[0] == "stale"


@pytest.mark.short
def test_validation_index_shows_worst_of_several_validators(tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    (out / "a.txt").write_text("x")
    _verdict(out, "a.txt", 0, i=0)
    _verdict(out, "a.txt", 1, i=1)
    _verdict(out, "a.txt", 127, i=2)
    assert verdict_for(out / "a.txt", load_verdicts(out))[0] == "fail"


@pytest.mark.short
def test_validation_index_stale_after_script_edit(tmp_path):
    from omnibenchmark.core.validators import script_sha256

    out = tmp_path / "out"
    out.mkdir()
    (out / "a.txt").write_text("x")
    script = tmp_path / "v.py"
    script.write_text("pass\n")
    result = out / ".validation" / "a.txt.0.json"
    result.parent.mkdir()
    result.write_text(
        json.dumps(
            {
                "output_file": "a.txt",
                "exit_code": 0,
                "script": str(script),
                "script_sha256": script_sha256(str(script)),
            }
        )
    )
    assert verdict_for(out / "a.txt", load_verdicts(out))[0] == "pass"
    script.write_text("raise SystemExit(1)\n")
    assert verdict_for(out / "a.txt", load_verdicts(out))[0] == "stale"
