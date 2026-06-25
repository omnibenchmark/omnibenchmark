"""Unit test for the describe report's validation read-side.

The validators run records results under the *human-readable symlink* path
(snakemake's glob skips the dot-prefixed hash dir), while the status report's
filedict references outputs by their *hash* path. Both resolve to the same real
file, so `_load_validation_index` / `_validation_for` must bridge the two.
"""

import json

import pytest

from omnibenchmark.cli.describe import _load_validation_index, _validation_for


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

    index = _load_validation_index(out)

    # both aliases of the same real file resolve to the recorded verdict
    assert _validation_for(real_file, index)[0] == "pass"
    assert _validation_for(symlink_dir / "D2_data.json", index)[0] == "pass"
    # the log path is carried through (absolute)
    assert _validation_for(real_file, index)[1].endswith("D2_data.json.log")
    # an output with no result -> not_validated
    assert _validation_for(out / "data" / "D1" / "x.json", index)[0] == "not_validated"
