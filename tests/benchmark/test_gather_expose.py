"""`gather[].expose`: one member per group passed as its own flag (#389).

The PCA case: every embedding is gathered per dataset, and the one whose
lineage carries `role: reference` is also passed as `--reference`.
"""

from types import SimpleNamespace

import pytest

from omnibenchmark.core._expand import expand_gather_stage
from omnibenchmark.model.benchmark import Benchmark
from omnibenchmark.model.resolved import TemplateContext


def _member(node_id, parent_id, stage_id, module_id, **labels):
    return SimpleNamespace(
        id=node_id,
        parent_id=parent_id,
        stage_id=stage_id,
        module_id=module_id,
        template_context=TemplateContext(provides={"name": module_id, **labels}),
    )


def _expand(pca_modules, expose, path_exclusions=None):
    """Two datasets x `pca_modules` [(module id, role)], gathered by data."""
    nodes_by_id = {}
    producers = []
    for d in ("d1", "d2"):
        nodes_by_id[d] = _member(d, None, "data", d)
        for mod, role in pca_modules:
            nid = f"{d}-pca-{mod}"
            nodes_by_id[nid] = _member(nid, d, "pca", mod, role=role)
            producers.append((nid, f"{d}/{mod}.emb"))
    stage = SimpleNamespace(
        id="agree",
        gather=[SimpleNamespace(from_="pca.embedding", group_by="data", expose=expose)],
        modules=[
            SimpleNamespace(
                id="subspace", name=None, parameters=None, provides=None, resources=None
            )
        ],
        outputs=[SimpleNamespace(id="agree.tsv", path="agree.tsv")],
        resources=None,
    )
    model = SimpleNamespace(
        get_name=lambda: "b", get_version=lambda: "1", get_author=lambda: "me"
    )
    return expand_gather_stage(
        stage=stage,
        benchmark=SimpleNamespace(model=model),
        resolved_modules_cache={("agree", "subspace"): object()},
        output_to_nodes={"pca.embedding": producers},
        nodes_by_id=nodes_by_id,
        path_exclusions=path_exclusions,
    )


def _flags(node):
    out = {}
    for key, flag in node.input_name_mapping.items():
        out.setdefault(flag, []).append(node.inputs[key])
    return out


PCA = [("sklearn", "reference"), ("irlba", "irlba"), ("rand", "rand")]


@pytest.mark.short
def test_expose_passes_the_reference_as_its_own_flag():
    nodes = _expand(PCA, {"reference": {"role": "reference"}})
    by_id = {n.id: n for n in nodes}
    d1 = _flags(by_id["agree-subspace-d1.default"])
    assert d1["reference"] == ["d1/sklearn.emb"]
    # Not a filter: the reference stays in the gathered list.
    assert d1["pca.embedding"] == ["d1/sklearn.emb", "d1/irlba.emb", "d1/rand.emb"]
    assert _flags(by_id["agree-subspace-d2.default"])["reference"] == ["d2/sklearn.emb"]
    # Provenance lists each member once.
    assert len(by_id["agree-subspace-d1.default"].gathered_from) == 3


@pytest.mark.short
def test_expose_by_builtin_name():
    nodes = _expand(PCA, {"reference": {"name": "irlba"}})
    assert _flags(nodes[0])["reference"] == ["d1/irlba.emb"]


@pytest.mark.short
def test_expose_without_a_match_fails_at_plan_time():
    with pytest.raises(ValueError, match="exactly one member.*found 0"):
        _expand(PCA, {"reference": {"role": "nope"}})


@pytest.mark.short
def test_expose_with_two_matches_fails_at_plan_time():
    two = PCA + [("exact", "reference")]
    with pytest.raises(ValueError, match="found 2: d1-pca-sklearn, d1-pca-exact"):
        _expand(two, {"reference": {"role": "reference"}})


@pytest.mark.short
def test_expose_counts_after_exclusion():
    # Excluding the reference for this gather module leaves nothing to expose.
    with pytest.raises(ValueError, match="found 0"):
        _expand(
            PCA,
            {"reference": {"role": "reference"}},
            path_exclusions={"subspace": ["sklearn"]},
        )


_YAML = """
id: t
description: t
version: '1.0'
benchmarker: me
api_version: {api}
software_backend: host
software_environments:
  env: {{description: e, easyconfig: e.eb}}
stages:
  - id: data
    outputs: [{{id: data.raw, path: d.json}}]
    modules: [{{id: d1, {repo}}}]
  - id: pca
    provides: [role]
    inputs: [data.raw]
    outputs: [{{id: pca.embedding, path: e.json}}]
    modules:
      - {{id: sklearn, provides: {{role: reference}}, {repo}}}
      - {{id: irlba, {repo}}}
  - id: agree
    gather:
      - {{from: pca.embedding, group_by: data, expose: {expose}}}
    outputs: [{{id: agree.tsv, path: agree.tsv}}]
    modules: [{{id: subspace, {repo}}}]
"""
_REPO = "repository: {url: 'http://x', commit: abc}, software_environment: env"


def _parse(api="0.8.0", expose="{reference: {role: reference}}"):
    return Benchmark.from_yaml(_YAML.format(api=api, expose=expose, repo=_REPO))


@pytest.mark.short
def test_expose_parses():
    spec = _parse().stages[2].gather[0]
    assert spec.expose == {"reference": {"role": "reference"}}


@pytest.mark.parametrize(
    "api, expose, message",
    [
        ("0.7.0", "{reference: {role: reference}}", "requires api_version ≥ 0.8.0"),
        ("0.8.0", "{pca.embedding: {role: reference}}", "clash"),
        ("0.8.0", "{reference: {rol: reference}}", "no stage declares"),
    ],
)
@pytest.mark.short
def test_expose_rejected_at_parse_time(api, expose, message):
    with pytest.raises(Exception, match=message):
        _parse(api=api, expose=expose)
