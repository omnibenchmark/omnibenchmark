"""Regression tests for `ob describe status ... --logs` (issue #318).

The status reporter used to look for ``stdout.log``/``stderr.log`` inside each
node's *output* directory, but the Snakemake backend actually writes a single
combined log to ``{out_dir}/.logs/{rule_name}.log`` where ``rule_name`` is the
sanitized, fully-nested node id.  As a result ``--logs`` never showed anything.

These tests pin the status side to the backend's naming and exercise the
log-selection logic in ``print_exec_path_dict``.
"""

from types import SimpleNamespace

from omnibenchmark.backend.snakemake import SnakemakeGenerator
from omnibenchmark.core._paths import sanitize_rule_name
from omnibenchmark.core.status.status import (
    _attach_log_paths,
    print_exec_path_dict,
)


def _node(stage_id, module_id, param_id):
    return SimpleNamespace(
        stage_id=stage_id,
        module_id=module_id,
        param_id=param_id,
        parameters=None,
    )


class _FakeStage:
    """Minimal stand-in for ExecutionPathStage."""

    def __init__(self, node, code, output_file="out/x"):
        self.node = node
        self._code = code
        self.log = None
        self._output_file = output_file

    def dict_encode(self, full=False, cumulative=False):
        return self._code

    def get_output_files(self, *a, **k):
        return [self._output_file]


class _FakeExecPath:
    """Minimal stand-in for ExecutionPath."""

    def __init__(self, stage_nodes):
        # stage_nodes: list of (stage_id, _FakeStage)
        self.stages = [s for s, _ in stage_nodes]
        self.exec_path = {s: st for s, st in stage_nodes}

    def dict_encode(self, full=False, cumulative=False, stages=None):
        codes = "".join(self.exec_path[s].dict_encode() for s in self.stages)
        nmiss = sum(1 for c in codes if c != ".")
        first_failed = next(
            (s for s in self.stages if self.exec_path[s].dict_encode() != "."),
            None,
        )
        return codes, nmiss, first_failed


def _build_path(specs):
    """specs: list of (stage, module, param_id, code) -> _FakeExecPath."""
    return _FakeExecPath(
        [(st, _FakeStage(_node(st, mod, pid), code)) for st, mod, pid, code in specs]
    )


def test_attach_log_paths_matches_backend_naming(tmp_path):
    """The reconstructed log file equals the backend's sanitized nested id."""
    logs_dir = tmp_path / ".logs"
    logs_dir.mkdir()

    # Two-stage path: a default-param dataset feeding a hashed-param method.
    nested_id = "data-D1.default-methods-M1.17e49005"
    rule_name = sanitize_rule_name(nested_id)
    # Sanity: a hashed param must collapse to a single underscore, not two.
    assert rule_name == "data_D1_default_methods_M1_17e49005"
    # Pin to the backend's own sanitizer (the source of the on-disk filename).
    assert rule_name == SnakemakeGenerator._sanitize_rule_name(
        SnakemakeGenerator.__new__(SnakemakeGenerator), nested_id
    )

    (logs_dir / "data_D1_default.log").write_text("data log")
    (logs_dir / f"{rule_name}.log").write_text("methods log")

    ep = _build_path(
        [("data", "D1", "default", "M"), ("methods", "M1", ".17e49005", "M")]
    )
    _attach_log_paths({0: ep}, tmp_path)

    # The intermediate stage uses its own (un-nested) id.
    assert ep.exec_path["data"].log == logs_dir / "data_D1_default.log"
    # The downstream stage uses the full nested id (single-underscore hash).
    assert ep.exec_path["methods"].log == logs_dir / f"{rule_name}.log"


def test_attach_log_paths_missing_file_is_none(tmp_path):
    (tmp_path / ".logs").mkdir()
    ep = _build_path([("data", "D1", "default", "M")])
    _attach_log_paths({0: ep}, tmp_path)
    # No log file written -> attribute stays None.
    assert ep.exec_path["data"].log is None


def test_print_logs_shows_every_missing_stage(tmp_path):
    """--logs lists a LOG line for *each* missing-output stage, not just one."""
    stages = ["data", "methods", "metrics"]
    ep = _build_path(
        [
            ("data", "D1", "default", "."),  # complete, no log expected
            ("methods", "M1", "default", "M"),  # missing
            ("metrics", "m1", "default", "M"),  # missing
        ]
    )
    # Give the two failing stages real log paths.
    ep.exec_path["methods"].log = tmp_path / "methods.log"
    ep.exec_path["metrics"].log = tmp_path / "metrics.log"

    out = print_exec_path_dict({0: ep}, stages, threshold_n_missing=1, logs=True)

    assert "LOG (methods):" in out
    assert "LOG (metrics):" in out
    # The complete stage has no log line.
    assert "LOG (data):" not in out


def test_print_logs_off_emits_nothing(tmp_path):
    stages = ["data", "methods"]
    ep = _build_path(
        [("data", "D1", "default", "."), ("methods", "M1", "default", "M")]
    )
    ep.exec_path["methods"].log = tmp_path / "methods.log"

    out = print_exec_path_dict({0: ep}, stages, threshold_n_missing=1, logs=False)
    assert "LOG (" not in out
