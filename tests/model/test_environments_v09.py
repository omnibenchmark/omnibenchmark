"""api 0.9.0 software environments: structured backends and an included file.

software_environments:            # inline ...
  py:
    - {type: conda, uri: envs/py.yaml}
    - {type: apptainer, uri: oras://registry/py.sif}

software_environments: environments.yaml   # ... or a path to that mapping
"""

import textwrap

import pytest
import yaml

from omnibenchmark.model import Benchmark
from omnibenchmark.model.validation import BenchmarkParseError

PY = [
    {"type": "conda", "uri": "envs/py.yaml"},
    {"type": "apptainer", "uri": "oras://registry/py.sif"},
]


def _doc(envs, api="0.9.0", backend="conda"):
    return {
        "id": "b",
        "benchmarker": "me",
        "version": "1.0",
        "api_version": api,
        "software_backend": backend,
        "software_environments": envs,
        "stages": [
            {
                "id": "data",
                "modules": [
                    {
                        "id": "D1",
                        "repository": {
                            "url": "https://example.org/d.git",
                            "commit": "abc1234",
                        },
                    }
                ],
                "outputs": [{"id": "data_counts", "path": "{dataset}.txt"}],
            }
        ],
    }


def _write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(data, sort_keys=False))
    return path


def _env(b, env_id="py"):
    return b.get_software_environments()[env_id]


@pytest.mark.short
class TestStructuredBackends:
    def test_list_maps_onto_backend_fields(self, tmp_path):
        b = Benchmark.from_yaml(_write(tmp_path / "plan.yaml", _doc({"py": PY})))
        env = _env(b)
        assert env.conda == "envs/py.yaml"
        assert env.apptainer == "oras://registry/py.sif"
        assert env.docker is None

    def test_envmodules_and_docker(self, tmp_path):
        envs = {
            "py": [
                {"type": "envmodules", "uri": "Python/3.12"},
                {"type": "docker", "uri": "python:3.12"},
            ]
        }
        env = _env(Benchmark.from_yaml(_write(tmp_path / "p.yaml", _doc(envs))))
        assert env.envmodule == "Python/3.12"
        assert env.docker == "python:3.12"

    def test_works_from_a_yaml_string(self):
        b = Benchmark.from_yaml(yaml.safe_dump(_doc({"py": PY})))
        assert _env(b).conda == "envs/py.yaml"

    @pytest.mark.parametrize(
        "entry, match",
        [
            ({"type": "pip", "uri": "x"}, "unknown backend type 'pip'"),
            ({"type": "host", "uri": "x"}, "unknown backend type 'host'"),
            ({"type": "conda"}, "needs `type` and `uri`"),
            ({"uri": "x"}, "needs `type` and `uri`"),
            ({"type": "conda", "uri": "x", "extra": 1}, "unexpected key"),
            ("conda", "needs `type` and `uri`"),
        ],
    )
    def test_bad_entry(self, entry, match):
        with pytest.raises(BenchmarkParseError, match=match):
            Benchmark.from_yaml(yaml.safe_dump(_doc({"py": [entry]})))

    def test_duplicate_type(self):
        envs = {"py": [PY[0], {"type": "conda", "uri": "envs/other.yaml"}]}
        with pytest.raises(
            BenchmarkParseError, match="'py' lists backend 'conda' twice"
        ):
            Benchmark.from_yaml(yaml.safe_dump(_doc(envs)))

    def test_legacy_mapping_rejected_at_0_9(self):
        envs = {"py": {"description": "x", "conda": "envs/py.yaml"}}
        with pytest.raises(BenchmarkParseError, match=r"list of \{type, uri\}"):
            Benchmark.from_yaml(yaml.safe_dump(_doc(envs)))

    def test_list_form_rejected_before_0_9(self):
        with pytest.raises(BenchmarkParseError, match="api_version ≥ 0.9.0"):
            Benchmark.from_yaml(yaml.safe_dump(_doc({"py": PY}, api="0.7.0")))

    def test_legacy_mapping_still_fine_before_0_9(self):
        envs = {"py": {"description": "x", "conda": "envs/py.yaml"}}
        b = Benchmark.from_yaml(yaml.safe_dump(_doc(envs, api="0.7.0")))
        assert _env(b).conda == "envs/py.yaml"

    def test_sole_env_auto_fill_still_applies(self):
        b = Benchmark.from_yaml(yaml.safe_dump(_doc({"py": PY})))
        assert b.stages[0].modules[0].software_environment == "py"


@pytest.mark.short
class TestIncludedFile:
    def test_include_next_to_plan(self, tmp_path):
        _write(tmp_path / "environments.yaml", {"py": PY})
        plan = _write(tmp_path / "plan.yaml", _doc("environments.yaml"))
        assert _env(Benchmark.from_yaml(plan)).conda == "envs/py.yaml"

    def test_local_uris_are_relative_to_the_included_file(self, tmp_path):
        _write(tmp_path / "shared" / "environments.yaml", {"py": PY})
        plan = _write(tmp_path / "plan.yaml", _doc("shared/environments.yaml"))
        env = _env(Benchmark.from_yaml(plan))
        assert env.conda == "shared/envs/py.yaml"
        assert env.apptainer == "oras://registry/py.sif"  # URLs untouched

    def test_absolute_uri_untouched(self, tmp_path):
        envs = {"py": [{"type": "conda", "uri": "/opt/envs/py.yaml"}]}
        _write(tmp_path / "shared" / "environments.yaml", envs)
        plan = _write(tmp_path / "plan.yaml", _doc("shared/environments.yaml"))
        assert _env(Benchmark.from_yaml(plan)).conda == "/opt/envs/py.yaml"

    def test_non_path_backends_untouched(self, tmp_path):
        envs = {
            "py": [
                {"type": "envmodules", "uri": "Python/3.12"},
                {"type": "docker", "uri": "python:3.12"},
            ]
        }
        _write(tmp_path / "shared" / "environments.yaml", envs)
        plan = _write(tmp_path / "plan.yaml", _doc("shared/environments.yaml"))
        env = _env(Benchmark.from_yaml(plan))
        assert env.envmodule == "Python/3.12"
        assert env.docker == "python:3.12"

    def test_hash_follows_included_contents(self, tmp_path):
        inc = _write(tmp_path / "environments.yaml", {"py": PY})
        plan = _write(tmp_path / "plan.yaml", _doc("environments.yaml"))
        before = Benchmark.from_yaml(plan).summary_hash()
        _write(inc, {"py": [{"type": "conda", "uri": "envs/py2.yaml"}]})
        assert Benchmark.from_yaml(plan).summary_hash() != before

    def test_hash_same_as_inline(self, tmp_path):
        _write(tmp_path / "environments.yaml", {"py": PY})
        included = Benchmark.from_yaml(
            _write(tmp_path / "a.yaml", _doc("environments.yaml"))
        )
        inline = Benchmark.from_yaml(_write(tmp_path / "b.yaml", _doc({"py": PY})))
        assert included.summary_hash() == inline.summary_hash()

    def test_missing_file(self, tmp_path):
        plan = _write(tmp_path / "plan.yaml", _doc("nope.yaml"))
        with pytest.raises(BenchmarkParseError, match="not found"):
            Benchmark.from_yaml(plan)

    def test_included_file_must_be_a_mapping(self, tmp_path):
        (tmp_path / "environments.yaml").write_text("- py\n")
        plan = _write(tmp_path / "plan.yaml", _doc("environments.yaml"))
        with pytest.raises(BenchmarkParseError, match="mapping of environment ids"):
            Benchmark.from_yaml(plan)

    def test_error_points_into_included_file(self, tmp_path):
        inc = tmp_path / "environments.yaml"
        inc.write_text(
            textwrap.dedent(
                """\
                py:
                  - type: conda
                    uri: envs/py.yaml
                  - type: pip
                    uri: x
                """
            )
        )
        plan = _write(tmp_path / "plan.yaml", _doc("environments.yaml"))
        with pytest.raises(BenchmarkParseError) as exc:
            Benchmark.from_yaml(plan)
        assert exc.value.yaml_file == inc
        assert exc.value.line_number == 4

    def test_include_needs_a_plan_file(self):
        with pytest.raises(BenchmarkParseError, match="loaded from a file"):
            Benchmark.from_yaml(yaml.safe_dump(_doc("environments.yaml")))

    def test_include_rejected_before_0_9(self, tmp_path):
        _write(tmp_path / "environments.yaml", {"py": PY})
        plan = _write(tmp_path / "plan.yaml", _doc("environments.yaml", api="0.7.0"))
        with pytest.raises(BenchmarkParseError, match="api_version ≥ 0.9.0"):
            Benchmark.from_yaml(plan)
