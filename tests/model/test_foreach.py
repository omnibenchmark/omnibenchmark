"""Model rules for `foreach` outputs (design 013 §3)."""

import pytest

from omnibenchmark.model.benchmark import Benchmark
from omnibenchmark.model.resolved import TemplateContext
from omnibenchmark.model.validation import BenchmarkParseError

_PLAN = """
id: b
description: d
version: "1.0"
benchmarker: t
api_version: {api}
software_backend: host
software_environments:
  - id: host
    description: host
stages:
  - id: data
    modules:
      - id: D1
        software_environment: host
        repository:
          url: https://github.com/t/t.git
          commit: abc123
    outputs:
{outputs}
"""

_OK = """      - id: counts
        path: "{sample}/counts.h5"
        foreach: sample"""


def _load(outputs=_OK, api="0.8.0"):
    return Benchmark.from_yaml(_PLAN.format(api=api, outputs=outputs))


@pytest.mark.short
def test_foreach_parses():
    stage = _load().stages[0]
    assert stage.outputs[0].foreach == "sample"


@pytest.mark.short
@pytest.mark.parametrize(
    "outputs, api, match",
    [
        (_OK, "0.7.0", "requires api_version"),
        (
            _OK.replace('"{sample}/counts.h5"', '"counts_{sample}.h5"'),
            "0.8.0",
            "must start",
        ),
        (
            _OK.replace('"{sample}/counts.h5"', '"{sample}/{sample}.h5"'),
            "0.8.0",
            "must start",
        ),
        (_OK.replace("foreach: sample", "foreach: params"), "0.8.0", "identifier"),
        (_OK + "\n" + _OK.replace("counts", "other"), "0.8.0", "at most one"),
    ],
)
def test_foreach_rejected(outputs, api, match):
    with pytest.raises(BenchmarkParseError, match=match):
        _load(outputs, api)


@pytest.mark.short
def test_substitute_keeps_capture():
    ctx = TemplateContext(provides={"sample": "nope"}, module_attrs={"id": "D1"})
    got = ctx.substitute("{sample}/{module.id}.h5", keep=("sample",))
    assert got == "{sample}/D1.h5"
    with pytest.raises(ValueError, match="Unresolved"):
        ctx.substitute("{other}/x.h5", keep=("sample",))


_OVER_PLAN = """
id: b
description: d
version: "1.0"
benchmarker: t
api_version: {api}
software_backend: host
software_environments:
  - id: host
    description: host
stages:
  - id: data
    modules:
      - id: D1
        software_environment: host
        repository: {{url: https://github.com/t/t.git, commit: abc123}}
    outputs:
      - id: counts
        path: "{{sample}}/counts.h5"
        foreach: sample
  - id: concat
    modules:
      - id: C
        software_environment: host
        repository: {{url: https://github.com/t/t.git, commit: abc123}}
{body}
    outputs:
      - id: all
        path: "all.h5"
"""


def _over(body, api="0.8.0"):
    return Benchmark.from_yaml(_OVER_PLAN.format(api=api, body=body))


@pytest.mark.short
def test_over_accepted():
    b = _over("    inputs: [counts]\n    over: sample")
    assert b.get_stage("concat").over == "sample"


@pytest.mark.short
@pytest.mark.parametrize(
    "body, api, match",
    [
        ("    inputs: [counts]\n    over: sample", "0.7.0", "requires api_version"),
        ("    inputs: [counts]\n    over: lane", "0.8.0", "names no `foreach` capture"),
        ("    over: sample", "0.8.0", "no inputs"),
        ("    inputs: [counts]\n    over: 'a-b'", "0.8.0", "must be an identifier"),
        (
            "    inputs: [counts]\n    over: sample\n"
            "    gather: [{from: counts, group_by: data}]",
            "0.8.0",
            "both `over` and `gather`",
        ),
    ],
)
def test_over_rejected(body, api, match):
    with pytest.raises(BenchmarkParseError, match=match):
        _over(body, api=api)


@pytest.mark.short
def test_over_absent_keeps_summary_hash():
    """`over: null` and no `over` hash the same: existing plans keep theirs."""
    a = _over("    inputs: [counts]").summary_hash()
    b = _over("    inputs: [counts]\n    over: null").summary_hash()
    assert a == b
