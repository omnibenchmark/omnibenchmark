"""Model rules for `foreach` outputs (design 012 §3)."""

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
