"""Regression tests for stale-repo encoding and cumulative invalid counts in status."""

import os

import pytest

from omnibenchmark.core.status.status import prepare_status, print_exec_path_dict
from omnibenchmark.model.repo import get_repo_hash

_TWO_STAGE_PLAN = """
id: stale-status-test
description: two-stage plan for stale status checks
version: "1.0"
benchmarker: test
api_version: "0.8.0"
software_backend: host
software_environments:
  host:
    description: host
stages:
  - id: data
    modules:
      - id: D1
        software_environment: host
        repository:
          url: https://example.org/data.git
          commit: aaa1111
    outputs:
      - id: data.counts
        path: "counts.txt"
  - id: methods
    modules:
      - id: M1
        software_environment: host
        repository:
          url: https://example.org/method.git
          commit: bbb2222
    inputs:
      - data.counts
    outputs:
      - id: methods.result
        path: "result.txt"
"""


def _write_with_mtime(path, content, mtime):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)
    os.utime(path, (mtime, mtime))


def _setup_stale_upstream_repo(tmp_path):
    """Create a 2-stage benchmark whose upstream repo is newer than its output."""
    plan = tmp_path / "plan.yaml"
    plan.write_text(_TWO_STAGE_PLAN)

    out_dir = tmp_path / "out"
    cache_dir = tmp_path / "repos"

    data_out = out_dir / "data.D1" / "counts.txt"
    method_out = out_dir / "data.D1" / "methods.M1" / "result.txt"
    data_repo_file = (
        cache_dir / get_repo_hash("https://example.org/data.git", "aaa1111") / "run.py"
    )
    method_repo_file = (
        cache_dir
        / get_repo_hash("https://example.org/method.git", "bbb2222")
        / "run.py"
    )

    # method output is newer than data output and method repo,
    # but data's repo was modified after data's output.
    _write_with_mtime(method_repo_file, "# method", 1000.0)
    _write_with_mtime(data_out, "counts", 2000.0)
    _write_with_mtime(method_out, "result", 3000.0)
    _write_with_mtime(data_repo_file, "# updated data repo", 4000.0)

    return plan, out_dir, cache_dir


@pytest.mark.short
def test_stale_repo_encodes_as_r_not_question_mark(tmp_path):
    """`dict_encode(full=True)` returns 'R' when a repo is newer than outputs."""
    plan, out_dir, cache_dir = _setup_stale_upstream_repo(tmp_path)

    _, _, exec_path_dict = prepare_status(
        str(plan), str(out_dir), return_all=True, cache_dir=cache_dir
    )

    code, nmiss, first_failed = exec_path_dict[0].dict_encode(
        full=True, cumulative=True
    )
    assert code == "RR"
    assert nmiss == 2
    assert first_failed == "data"

    rendered = print_exec_path_dict(
        exec_path_dict, ["data", "methods"], threshold_n_missing=1, full=True
    )
    assert "RR" in rendered
    assert "?" not in rendered


@pytest.mark.short
def test_cumulative_filedict_counts_downstream_invalid_outputs(tmp_path):
    """`create_filedict(cumulative=True)` propagates invalid counts downstream."""
    plan, out_dir, cache_dir = _setup_stale_upstream_repo(tmp_path)

    status_dict, filedict, _ = prepare_status(
        str(plan), str(out_dir), return_all=True, cache_dir=cache_dir
    )

    assert status_dict["stages"]["data"]["n_invalid"] == 1
    assert status_dict["stages"]["data"]["n_invalid_repo_is_newer"] == 1
    assert status_dict["stages"]["methods"]["n_observed"] == 0
    assert status_dict["stages"]["methods"]["n_invalid_repo_is_newer"] == 1
    assert status_dict["stages"]["methods"]["n_invalid"] == 1
    (method_node,) = filedict["methods"].keys()
    assert filedict["methods"][method_node]["n_invalid"] == 1


@pytest.mark.short
def test_stale_input_file_reports_f_code(tmp_path):
    """An input file newer than downstream outputs encodes as 'F'."""
    plan = tmp_path / "plan.yaml"
    plan.write_text(_TWO_STAGE_PLAN)

    out_dir = tmp_path / "out"
    cache_dir = tmp_path / "repos"

    data_out = out_dir / "data.D1" / "counts.txt"
    method_out = out_dir / "data.D1" / "methods.M1" / "result.txt"
    data_repo_file = (
        cache_dir / get_repo_hash("https://example.org/data.git", "aaa1111") / "run.py"
    )
    method_repo_file = (
        cache_dir
        / get_repo_hash("https://example.org/method.git", "bbb2222")
        / "run.py"
    )

    _write_with_mtime(data_repo_file, "# data", 1000.0)
    _write_with_mtime(method_repo_file, "# method", 1000.0)
    _write_with_mtime(method_out, "stale result", 2000.0)
    _write_with_mtime(data_out, "fresh counts", 3000.0)

    status_dict, _, exec_path_dict = prepare_status(
        str(plan), str(out_dir), return_all=True, cache_dir=cache_dir
    )

    code, nmiss, first_failed = exec_path_dict[0].dict_encode(
        full=True, cumulative=True
    )
    assert code == ".F"
    assert nmiss == 1
    assert first_failed == "methods"
    assert status_dict["stages"]["methods"]["n_invalid"] == 1
    assert status_dict["stages"]["methods"]["n_invalid_input_file_is_newer"] == 1
