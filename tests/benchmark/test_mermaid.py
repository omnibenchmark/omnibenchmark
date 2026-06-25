"""Topology-focused tests for Mermaid diagram generation.

These build benchmark *models* directly (no execution DAG / out_dir needed) and
assert on the module-to-module edges the generator emits. They cover the
regression where an input group spanning multiple upstream stages dropped every
edge except the one to the topologically-latest stage.
"""

import re

from omnibenchmark.model import Benchmark
from omnibenchmark.core._mermaid import generate_mermaid_diagram

_HEADER = """\
id: {id}
description: stress
version: "1.0"
benchmarker: x
storage: https://storage.github.com/
api_version: "0.3.0"
storage_api: S3
storage_bucket_name: b
software_backend: host
software_environments:
  e:
    description: e
    easyconfig: e.eb
    envmodule: e
    conda: e.yaml
    apptainer: http://registry.ch/x.sif
"""

_REPO = "{url: https://github.com/omnibenchmark-example/data.git, commit: 63b7b36}"


def _module(mid, params=None):
    """*params* is a list of parameter sets. A set given as a list renders the
    deprecated ``values:`` form; given as a dict it renders the dict form (the
    only form ``--compact-params`` can display)."""
    lines = [
        f"      - id: {mid}",
        "        software_environment: e",
        f"        repository: {_REPO}",
    ]
    if params:
        lines.append("        parameters:")
        for p in params:
            if isinstance(p, dict):
                first = True
                for k, v in p.items():
                    prefix = "          - " if first else "            "
                    lines.append(f"{prefix}{k}: {v}")
                    first = False
            else:
                vals = ", ".join(f'"{v}"' for v in p)
                lines.append(f"          - values: [{vals}]")
    return "\n".join(lines)


def _stage(sid, modules, inputs=None, output_id=None):
    """Compose a stage block. *inputs* is a list of input groups; each group is
    a list of output ids (flat-list form). *modules* is a list of (id, params)."""
    out = output_id or f"{sid}.out"
    block = [f"  - id: {sid}", "    modules:"]
    block += [_module(mid, params) for mid, params in modules]
    if inputs is not None:
        block.append("    inputs:")
        for group in inputs:
            if len(group) == 1:
                block.append(f"      - {group[0]}")
            else:
                block.append("      - entries:")
                block += [f"          - {e}" for e in group]
    block.append("    outputs:")
    block.append(f'      - {{id: {out}, path: "{sid}.txt"}}')
    return "\n".join(block)


def _model(bid, stages):
    yaml = _HEADER.format(id=bid) + "stages:\n" + "\n".join(stages)
    return Benchmark.from_yaml(yaml)


def _edges(model, **kw):
    diagram = generate_mermaid_diagram(model, show_params=False, **kw)
    return {(a, b) for a, b in re.findall(r"(\w+) --> (\w+)", diagram)}


# --------------------------------------------------------------------------- #
# Edge topology
# --------------------------------------------------------------------------- #


def test_single_input_group_spanning_two_stages():
    """Regression: one input group that references outputs from two upstream
    stages must produce edges from BOTH, not only the latest one."""
    model = _model(
        "diamond",
        [
            _stage("data", [("D", None)], output_id="data.out"),
            _stage("proc", [("P", None)], inputs=[["data.out"]], output_id="proc.out"),
            _stage(
                "combine",
                [("C", None)],
                inputs=[["data.out", "proc.out"]],
                output_id="combine.out",
            ),
        ],
    )
    assert _edges(model) == {("D", "P"), ("D", "C"), ("P", "C")}


def test_input_group_spanning_three_stages():
    model = _model(
        "triple",
        [
            _stage("data", [("D", None)], output_id="data.out"),
            _stage("s1", [("S1", None)], inputs=[["data.out"]], output_id="s1.out"),
            _stage("s2", [("S2", None)], inputs=[["s1.out"]], output_id="s2.out"),
            _stage(
                "fin",
                [("F", None)],
                inputs=[["data.out", "s1.out", "s2.out"]],
                output_id="fin.out",
            ),
        ],
    )
    assert _edges(model) == {
        ("D", "S1"),
        ("S1", "S2"),
        ("D", "F"),
        ("S1", "F"),
        ("S2", "F"),
    }


def test_separate_single_stage_groups_equivalent():
    """The deprecated multi-`entries` form (separate groups) yields the same
    edges as a single multi-stage group."""
    model = _model(
        "groups",
        [
            _stage("data", [("D", None)], output_id="data.out"),
            _stage("proc", [("P", None)], inputs=[["data.out"]], output_id="proc.out"),
            _stage(
                "combine",
                [("C", None)],
                inputs=[["data.out"], ["proc.out"]],
                output_id="combine.out",
            ),
        ],
    )
    assert _edges(model) == {("D", "P"), ("D", "C"), ("P", "C")}


def test_skip_level_dependency():
    """A stage depending only on an early stage skips the intermediate one."""
    model = _model(
        "skip",
        [
            _stage("data", [("D", None)], output_id="data.out"),
            _stage("proc", [("P", None)], inputs=[["data.out"]], output_id="proc.out"),
            _stage("tail", [("T", None)], inputs=[["data.out"]], output_id="tail.out"),
        ],
    )
    assert _edges(model) == {("D", "P"), ("D", "T")}


def test_fan_in_from_independent_roots():
    """Two independent initial stages feeding one downstream stage both connect."""
    model = _model(
        "fanin",
        [
            _stage("a", [("A", None)], output_id="a.out"),
            _stage("b", [("B", None)], output_id="b.out"),
            _stage(
                "merge",
                [("M", None)],
                inputs=[["a.out", "b.out"]],
                output_id="merge.out",
            ),
        ],
    )
    assert _edges(model) == {("A", "M"), ("B", "M")}


def test_full_cartesian_across_modules():
    """Every upstream module connects to every downstream module in a stage."""
    model = _model(
        "cart",
        [
            _stage("data", [("D1", None), ("D2", None)], output_id="data.out"),
            _stage(
                "proc",
                [("P1", None), ("P2", None)],
                inputs=[["data.out"]],
                output_id="proc.out",
            ),
        ],
    )
    assert _edges(model) == {
        ("D1", "P1"),
        ("D1", "P2"),
        ("D2", "P1"),
        ("D2", "P2"),
    }


def test_initial_stage_has_no_incoming_edges():
    model = _model(
        "init",
        [
            _stage("data", [("D", None)], output_id="data.out"),
            _stage("proc", [("P", None)], inputs=[["data.out"]], output_id="proc.out"),
        ],
    )
    edges = _edges(model)
    assert not any(target == "D" for _, target in edges)


def test_no_self_edges():
    model = _model(
        "self",
        [
            _stage("data", [("D", None)], output_id="data.out"),
            _stage("proc", [("P", None)], inputs=[["data.out"]], output_id="proc.out"),
        ],
    )
    assert not any(a == b for a, b in _edges(model))


# --------------------------------------------------------------------------- #
# Structure / params
# --------------------------------------------------------------------------- #


def test_header_and_declaration():
    model = _model("MyBench", [_stage("data", [("D", None)], output_id="data.out")])
    diagram = generate_mermaid_diagram(model, show_params=False)
    assert "title: MyBench" in diagram
    assert "flowchart LR" in diagram
    assert "classDef param fill:#f96" in diagram


def test_stage_subgraphs_present():
    model = _model(
        "subg",
        [
            _stage("data", [("D", None)], output_id="data.out"),
            _stage("proc", [("P", None)], inputs=[["data.out"]], output_id="proc.out"),
        ],
    )
    diagram = generate_mermaid_diagram(model, show_params=False)
    assert "subgraph data" in diagram
    assert "subgraph proc" in diagram


def test_show_params_toggle():
    model = _model(
        "pbench",
        [
            _stage("data", [("D", None)], output_id="data.out"),
            _stage(
                "proc",
                [("P", [["-a 0", "-b 1"]])],
                inputs=[["data.out"]],
                output_id="proc.out",
            ),
        ],
    )
    assert "subgraph params_P" not in generate_mermaid_diagram(model, show_params=False)
    assert "subgraph params_P" in generate_mermaid_diagram(model, show_params=True)


def test_compact_params_render_one_node_per_paramset():
    model = _model(
        "cbench",
        [
            _stage("data", [("D", None)], output_id="data.out"),
            _stage(
                "proc",
                [("P", [{"alpha": 0.1}, {"alpha": 0.5}])],
                inputs=[["data.out"]],
                output_id="proc.out",
            ),
        ],
    )
    diagram = generate_mermaid_diagram(model, show_params=True, compact_params=True)
    assert "P_paramset0" in diagram
    assert "P_paramset1" in diagram
    # The param subgraph is wired to its module
    assert "params_P:::param --o P" in diagram


def test_module_without_params_has_no_param_subgraph():
    model = _model(
        "nopar",
        [
            _stage("data", [("D", None)], output_id="data.out"),
            _stage("proc", [("P", None)], inputs=[["data.out"]], output_id="proc.out"),
        ],
    )
    diagram = generate_mermaid_diagram(model, show_params=True)
    assert "subgraph params_" not in diagram


def test_graph_argument_is_ignored():
    """Passing an (irrelevant) graph object must not change the output."""
    model = _model(
        "ign",
        [
            _stage("data", [("D", None)], output_id="data.out"),
            _stage("proc", [("P", None)], inputs=[["data.out"]], output_id="proc.out"),
        ],
    )
    assert generate_mermaid_diagram(
        model, graph=None, show_params=False
    ) == generate_mermaid_diagram(model, graph=object(), show_params=False)
