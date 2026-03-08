"""Short unit tests for cli/run.py pure and near-pure helper functions."""

import pytest
from pathlib import Path
from unittest.mock import MagicMock, patch
from pydantic import ValidationError as PydanticValidationError

from click.testing import CliRunner
from omnibenchmark.cli.run import (
    format_pydantic_errors,
    _read_rule_log,
    _run_benchmark,
    _run_snakemake,
    _populate_git_cache,
    _substitute_params_in_path,
    _build_template_context,
    _select_input_nodes,
    _satisfies_requires,
    run,
)
from omnibenchmark.model import SoftwareBackendEnum
from omnibenchmark.model.resolved import TemplateContext
from omnibenchmark.model.params import Params


# ---------------------------------------------------------------------------
# format_pydantic_errors
# ---------------------------------------------------------------------------


def _make_pydantic_error(errors: list[dict]):
    """Create a PydanticValidationError from a list of raw error dicts."""
    from pydantic import BaseModel

    class _M(BaseModel):
        x: int

    try:
        _M(x="bad")
    except PydanticValidationError:
        # Patch internal errors for custom test data
        pass

    mock_err = MagicMock(spec=PydanticValidationError)
    mock_err.errors.return_value = errors
    return mock_err


@pytest.mark.short
class TestFormatPydanticErrors:
    def test_missing_field(self):
        err = _make_pydantic_error(
            [{"loc": ("stages",), "msg": "Field required", "type": "missing"}]
        )
        result = format_pydantic_errors(err)
        assert "Missing required field: 'stages'" in result

    def test_non_missing_field(self):
        err = _make_pydantic_error(
            [
                {
                    "loc": ("version",),
                    "msg": "Input should be a string",
                    "type": "string_type",
                }
            ]
        )
        result = format_pydantic_errors(err)
        assert "Field 'version'" in result
        assert "Input should be a string" in result

    def test_multiple_errors(self):
        err = _make_pydantic_error(
            [
                {"loc": ("id",), "msg": "Field required", "type": "missing"},
                {"loc": ("version",), "msg": "bad value", "type": "value_error"},
            ]
        )
        result = format_pydantic_errors(err)
        assert "Validation failed:" in result
        assert "Missing required field: 'id'" in result
        assert "Field 'version'" in result

    def test_nested_loc(self):
        err = _make_pydantic_error(
            [{"loc": ("stages", 0, "id"), "msg": "bad", "type": "value_error"}]
        )
        result = format_pydantic_errors(err)
        assert "stages -> 0 -> id" in result


# ---------------------------------------------------------------------------
# _read_rule_log
# ---------------------------------------------------------------------------


@pytest.mark.short
class TestReadRuleLog:
    def test_returns_none_when_missing(self, tmp_path):
        assert _read_rule_log(tmp_path, "some_rule") is None

    def test_reads_existing_log(self, tmp_path):
        log_dir = tmp_path / ".logs"
        log_dir.mkdir()
        (log_dir / "my_rule.log").write_text("some log content")
        result = _read_rule_log(tmp_path, "my_rule")
        assert result == "some log content"

    def test_returns_none_on_read_error(self, tmp_path):
        log_dir = tmp_path / ".logs"
        log_dir.mkdir()
        log_path = log_dir / "bad_rule.log"
        log_path.write_text("x")
        with patch.object(Path, "read_text", side_effect=OSError("perm")):
            result = _read_rule_log(tmp_path, "bad_rule")
        assert result is None


# ---------------------------------------------------------------------------
# _substitute_params_in_path
# ---------------------------------------------------------------------------


@pytest.mark.short
class TestSubstituteParamsInPath:
    def test_no_placeholder_unchanged(self):
        assert (
            _substitute_params_in_path("data/D1/out.json", None) == "data/D1/out.json"
        )

    def test_none_params_unchanged(self):
        assert (
            _substitute_params_in_path("data/{params.key}/out.json", None)
            == "data/{params.key}/out.json"
        )

    def test_substitutes_key(self):
        p = Params({"key": "value123"})
        result = _substitute_params_in_path("data/{params.key}/out.json", p)
        assert result == "data/value123/out.json"

    def test_missing_key_left_as_is(self):
        p = Params({"other": "x"})
        result = _substitute_params_in_path("data/{params.missing}/out.json", p)
        assert result == "data/{params.missing}/out.json"

    def test_multiple_substitutions(self):
        p = Params({"a": "1", "b": "2"})
        result = _substitute_params_in_path("{params.a}_{params.b}.txt", p)
        assert result == "1_2.txt"

    def test_no_params_placeholder_skips_regex(self):
        p = Params({"k": "v"})
        assert _substitute_params_in_path("plain/path.txt", p) == "plain/path.txt"


# ---------------------------------------------------------------------------
# _build_template_context
# ---------------------------------------------------------------------------


def _make_stage(stage_id, provides=None):
    s = MagicMock()
    s.id = stage_id
    s.provides = provides
    return s


def _make_input_node(module_id, stage_id, template_context=None):
    n = MagicMock()
    n.module_id = module_id
    n.stage_id = stage_id
    n.template_context = template_context
    return n


@pytest.mark.short
class TestBuildTemplateContext:
    def test_root_node_default_dataset(self):
        stage = _make_stage("data", provides=None)
        ctx = _build_template_context(stage, "D1")
        assert ctx.provides["dataset"] == "D1"
        assert ctx.module_attrs["id"] == "D1"
        assert ctx.module_attrs["stage"] == "data"

    def test_root_node_dataset_from_params(self):
        stage = _make_stage("data", provides=None)
        p = Params({"dataset": "pbmc3k"})
        ctx = _build_template_context(stage, "D1", params=p)
        assert ctx.provides["dataset"] == "pbmc3k"

    def test_root_node_provides_label_from_params(self):
        stage = _make_stage("data", provides=["treatment"])
        p = Params({"treatment": "ctrl"})
        ctx = _build_template_context(stage, "D1", params=p)
        assert ctx.provides["treatment"] == "ctrl"

    def test_root_node_provides_label_defaults_to_module_id(self):
        stage = _make_stage("data", provides=["treatment"])
        ctx = _build_template_context(stage, "D1")
        assert ctx.provides["treatment"] == "D1"

    def test_child_node_inherits_parent_context(self):
        parent_ctx = TemplateContext(
            provides={"dataset": "pbmc3k"},
            module_attrs={"id": "D1", "stage": "data"},
        )
        input_node = _make_input_node("D1", "data", template_context=parent_ctx)
        stage = _make_stage("methods", provides=None)
        ctx = _build_template_context(stage, "M1", input_node=input_node)
        assert ctx.provides["dataset"] == "pbmc3k"
        assert ctx.module_attrs["parent.id"] == "D1"
        assert ctx.module_attrs["parent.stage"] == "data"

    def test_child_node_no_parent_context(self):
        input_node = _make_input_node("D1", "data", template_context=None)
        stage = _make_stage("methods", provides=None)
        ctx = _build_template_context(stage, "M1", input_node=input_node)
        assert "dataset" not in ctx.provides
        assert ctx.module_attrs["parent.id"] == "D1"

    def test_child_node_adds_provides_label(self):
        parent_ctx = TemplateContext(
            provides={"dataset": "pbmc3k"},
            module_attrs={"id": "D1", "stage": "data"},
        )
        input_node = _make_input_node("D1", "data", template_context=parent_ctx)
        stage = _make_stage("methods", provides=["method"])
        ctx = _build_template_context(stage, "M1", input_node=input_node)
        assert ctx.provides["method"] == "M1"

    def test_child_node_provides_label_from_params(self):
        parent_ctx = TemplateContext(
            provides={"dataset": "pbmc3k"},
            module_attrs={"id": "D1", "stage": "data"},
        )
        input_node = _make_input_node("D1", "data", template_context=parent_ctx)
        stage = _make_stage("methods", provides=["method"])
        p = Params({"method": "kmeans"})
        ctx = _build_template_context(stage, "M1", input_node=input_node, params=p)
        assert ctx.provides["method"] == "kmeans"


# ---------------------------------------------------------------------------
# _satisfies_requires
# ---------------------------------------------------------------------------


@pytest.mark.short
class TestSatisfiesRequires:
    def test_no_template_context_returns_false(self):
        node = _make_input_node("D1", "data", template_context=None)
        assert _satisfies_requires({"dataset": "pbmc3k"}, node) is False

    def test_matching_single_constraint(self):
        ctx = TemplateContext(
            provides={"dataset": "pbmc3k"},
            module_attrs={"id": "D1", "stage": "data"},
        )
        node = _make_input_node("D1", "data", template_context=ctx)
        assert _satisfies_requires({"dataset": "pbmc3k"}, node) is True

    def test_mismatched_constraint(self):
        ctx = TemplateContext(
            provides={"dataset": "pbmc3k"},
            module_attrs={"id": "D1", "stage": "data"},
        )
        node = _make_input_node("D1", "data", template_context=ctx)
        assert _satisfies_requires({"dataset": "other"}, node) is False

    def test_missing_label_returns_false(self):
        ctx = TemplateContext(
            provides={},
            module_attrs={"id": "D1", "stage": "data"},
        )
        node = _make_input_node("D1", "data", template_context=ctx)
        assert _satisfies_requires({"dataset": "pbmc3k"}, node) is False

    def test_empty_requires_returns_true(self):
        ctx = TemplateContext(
            provides={"dataset": "pbmc3k"},
            module_attrs={"id": "D1", "stage": "data"},
        )
        node = _make_input_node("D1", "data", template_context=ctx)
        assert _satisfies_requires({}, node) is True

    def test_multiple_constraints_all_match(self):
        ctx = TemplateContext(
            provides={"dataset": "pbmc3k", "treatment": "ctrl"},
            module_attrs={"id": "D1", "stage": "data"},
        )
        node = _make_input_node("D1", "data", template_context=ctx)
        assert (
            _satisfies_requires({"dataset": "pbmc3k", "treatment": "ctrl"}, node)
            is True
        )

    def test_multiple_constraints_partial_mismatch(self):
        ctx = TemplateContext(
            provides={"dataset": "pbmc3k", "treatment": "ctrl"},
            module_attrs={"id": "D1", "stage": "data"},
        )
        node = _make_input_node("D1", "data", template_context=ctx)
        assert (
            _satisfies_requires({"dataset": "pbmc3k", "treatment": "stim"}, node)
            is False
        )


# ---------------------------------------------------------------------------
# _select_input_nodes
# ---------------------------------------------------------------------------


def _make_node(node_id, stage_id):
    n = MagicMock()
    n.id = node_id
    n.stage_id = stage_id
    return n


@pytest.mark.short
class TestSelectInputNodes:
    def test_empty_declared_inputs_returns_previous(self):
        prev = [_make_node("n1", "data")]
        result = _select_input_nodes([], {}, [], [], prev)
        assert result is prev

    def test_no_matching_outputs_returns_previous(self):
        prev = [_make_node("n1", "data")]
        result = _select_input_nodes(["data.raw"], {}, [], ["data"], prev)
        assert result is prev

    def test_single_input_selects_correct_stage(self):
        node_data = _make_node("data-D1", "data")
        node_methods = _make_node("methods-M1", "methods")
        resolved = [node_data, node_methods]
        output_to_nodes = {"data.raw": [("data-D1", "data/D1/out.json")]}
        stage_ids = ["data", "methods"]
        prev = []
        result = _select_input_nodes(
            ["data.raw"], output_to_nodes, resolved, stage_ids, prev
        )
        assert all(n.stage_id == "data" for n in result)

    def test_selects_deepest_providing_stage(self):
        n_data = _make_node("data-D1", "data")
        n_prep = _make_node("prep-P1", "preprocessing")
        resolved = [n_data, n_prep]
        output_to_nodes = {
            "data.raw": [("data-D1", "p1.json")],
            "prep.out": [("prep-P1", "p2.json")],
        }
        stage_ids = ["data", "preprocessing", "methods"]
        prev = []
        result = _select_input_nodes(
            ["data.raw", "prep.out"], output_to_nodes, resolved, stage_ids, prev
        )
        assert all(n.stage_id == "preprocessing" for n in result)

    def test_node_not_in_resolved_skipped(self):
        output_to_nodes = {"data.raw": [("ghost-node", "p.json")]}
        prev = [_make_node("fallback", "data")]
        result = _select_input_nodes(["data.raw"], output_to_nodes, [], ["data"], prev)
        assert result is prev


# ---------------------------------------------------------------------------
# _run_snakemake — missing Snakefile early exit
# ---------------------------------------------------------------------------


@pytest.mark.short
def test_run_snakemake_missing_snakefile_exits(tmp_path):
    """_run_snakemake should call sys.exit when Snakefile is absent."""

    with pytest.raises(SystemExit):
        _run_snakemake(
            out_dir=tmp_path,  # no Snakefile here
            cores=1,
            continue_on_error=False,
            software_backend=SoftwareBackendEnum.host,
            debug=False,
        )


# ---------------------------------------------------------------------------
# run command — dirty / unpinned warnings
# ---------------------------------------------------------------------------


@pytest.mark.short
def test_run_dirty_flag_logs_warning():
    with patch("omnibenchmark.cli.run._run_benchmark") as mock_rb:
        runner = CliRunner()
        runner.invoke(run, ["tests/data/mock_benchmark.yaml", "--dirty"])
        mock_rb.assert_called_once()
        assert mock_rb.call_args.kwargs.get("dirty") is True


@pytest.mark.short
def test_run_unpinned_flag_logs_warning():
    with patch("omnibenchmark.cli.run._run_benchmark") as mock_rb:
        runner = CliRunner()
        runner.invoke(run, ["tests/data/mock_benchmark.yaml", "--unpinned"])
        mock_rb.assert_called_once()
        assert mock_rb.call_args.kwargs.get("unpinned") is True


# ---------------------------------------------------------------------------
# _run_benchmark — PydanticValidationError + dry mode + remote storage
# ---------------------------------------------------------------------------


@pytest.mark.short
def test_run_benchmark_pydantic_error_exits(tmp_path):
    with patch("omnibenchmark.cli.run.BenchmarkExecution") as mock_be:
        mock_be.side_effect = PydanticValidationError.from_exception_data(
            title="Benchmark",
            line_errors=[],
        )
        with pytest.raises(SystemExit):
            _run_benchmark(
                benchmark_path="tests/data/mock_benchmark.yaml",
                cores=1,
                dry=False,
                continue_on_error=False,
                out_dir=str(tmp_path),
                debug=False,
                dirty=False,
            )


@pytest.mark.short
def test_run_benchmark_dry_host_exits_zero(tmp_path):
    with (
        patch("omnibenchmark.cli.run.BenchmarkExecution") as mock_be,
        patch("omnibenchmark.cli.run._populate_git_cache"),
        patch("omnibenchmark.cli.run._generate_explicit_snakefile"),
        patch("omnibenchmark.cli.run.write_run_manifest"),
    ):
        mock_b = MagicMock()
        mock_b.get_benchmark_software_backend.return_value = SoftwareBackendEnum.host
        mock_be.return_value = mock_b
        with pytest.raises(SystemExit) as exc:
            _run_benchmark(
                benchmark_path="tests/data/mock_benchmark.yaml",
                cores=1,
                dry=True,
                continue_on_error=False,
                out_dir=str(tmp_path),
                debug=False,
                dirty=False,
            )
        assert exc.value.code == 0


@pytest.mark.short
def test_run_benchmark_dry_conda_hint(tmp_path):
    with (
        patch("omnibenchmark.cli.run.BenchmarkExecution") as mock_be,
        patch("omnibenchmark.cli.run._populate_git_cache"),
        patch("omnibenchmark.cli.run._generate_explicit_snakefile"),
        patch("omnibenchmark.cli.run.write_run_manifest"),
        patch("omnibenchmark.cli.run.logger") as mock_logger,
    ):
        mock_b = MagicMock()
        mock_b.get_benchmark_software_backend.return_value = SoftwareBackendEnum.conda
        mock_be.return_value = mock_b
        with pytest.raises(SystemExit):
            _run_benchmark(
                benchmark_path="tests/data/mock_benchmark.yaml",
                cores=1,
                dry=True,
                continue_on_error=False,
                out_dir=str(tmp_path),
                debug=False,
                dirty=False,
            )
        # --use-conda hint should appear in the log
        logged = " ".join(str(c) for c in mock_logger.info.call_args_list)
        assert "--use-conda" in logged


@pytest.mark.short
def test_run_benchmark_remote_storage_builds_args(tmp_path):
    (tmp_path / "Snakefile").write_text("rule all: input: []")
    with (
        patch("omnibenchmark.cli.run.BenchmarkExecution") as mock_be,
        patch("omnibenchmark.cli.run._populate_git_cache"),
        patch("omnibenchmark.cli.run._generate_explicit_snakefile"),
        patch("omnibenchmark.cli.run.write_run_manifest"),
        patch("omnibenchmark.cli.run._run_snakemake") as mock_snakemake,
        patch("omnibenchmark.remote.storage.get_storage_from_benchmark"),
        patch(
            "omnibenchmark.remote.storage.remote_storage_snakemake_args"
        ) as mock_storage_args,
    ):
        mock_b = MagicMock()
        mock_b.get_benchmark_software_backend.return_value = SoftwareBackendEnum.host
        mock_be.return_value = mock_b
        mock_storage_args.return_value = {
            "default-storage-provider": "s3",
            "storage-s3-endpoint-url": "http://localhost:9000",
            "shared-fs-usage": None,  # None values should be skipped
            "use-singularity": False,  # False bool should be skipped
        }
        _run_benchmark(
            benchmark_path="tests/data/mock_benchmark.yaml",
            cores=1,
            dry=False,
            continue_on_error=False,
            out_dir=str(tmp_path),
            debug=False,
            dirty=False,
            use_remote_storage=True,
        )
        _, kwargs = mock_snakemake.call_args
        extra = kwargs["extra_snakemake_args"]
        assert "--default-storage-provider" in extra
        assert "s3" in extra
        assert "--storage-s3-endpoint-url" in extra
        # None and False values should NOT appear
        assert "--shared-fs-usage" not in extra
        assert "--use-singularity" not in extra


# ---------------------------------------------------------------------------
# _populate_git_cache — various paths
# ---------------------------------------------------------------------------


def _mock_benchmark(repo_url=None, commit="abc1234", local=False):
    """Build a minimal mock BenchmarkExecution for _populate_git_cache."""
    mock_module = MagicMock()
    if repo_url:
        mock_module.repository = MagicMock()
        mock_module.repository.url = repo_url
        mock_module.repository.commit = commit
    else:
        mock_module.repository = None

    mock_stage = MagicMock()
    mock_stage.modules = [mock_module]

    mock_b = MagicMock()
    mock_b.model.stages = [mock_stage]
    mock_b.model.metric_collectors = []
    return mock_b


@pytest.mark.short
def test_populate_git_cache_no_repos_quiet(tmp_path):
    """No repos → early return, no exception."""
    with patch("omnibenchmark.cli.run.get_git_cache_dir", return_value=tmp_path):
        _populate_git_cache(_mock_benchmark(), quiet=True, cores=1)


@pytest.mark.short
def test_populate_git_cache_no_repos_verbose(tmp_path):
    with patch("omnibenchmark.cli.run.get_git_cache_dir", return_value=tmp_path):
        _populate_git_cache(_mock_benchmark(), quiet=False, cores=1)


@pytest.mark.short
def test_populate_git_cache_local_path_skipped(tmp_path):
    """Local path repos are skipped → treated as no remote repos."""
    with (
        patch("omnibenchmark.cli.run.get_git_cache_dir", return_value=tmp_path),
        patch("omnibenchmark.cli.run.is_local_path", return_value=True),
    ):
        _populate_git_cache(
            _mock_benchmark(repo_url="/local/path/module"), quiet=True, cores=1
        )


@pytest.mark.short
def test_populate_git_cache_remote_repo_success(tmp_path):
    with (
        patch("omnibenchmark.cli.run.get_git_cache_dir", return_value=tmp_path),
        patch("omnibenchmark.cli.run.is_local_path", return_value=False),
        patch("omnibenchmark.cli.run.get_or_update_cached_repo") as mock_fetch,
        patch(
            "omnibenchmark.git.cache.parse_repo_url", return_value="github.com/org/repo"
        ),
    ):
        _populate_git_cache(
            _mock_benchmark(repo_url="https://github.com/org/repo.git"),
            quiet=False,
            cores=1,
        )
        mock_fetch.assert_called_once()


@pytest.mark.short
def test_populate_git_cache_remote_repo_failure(tmp_path):
    with (
        patch("omnibenchmark.cli.run.get_git_cache_dir", return_value=tmp_path),
        patch("omnibenchmark.cli.run.is_local_path", return_value=False),
        patch(
            "omnibenchmark.cli.run.get_or_update_cached_repo",
            side_effect=RuntimeError("network error"),
        ),
        patch(
            "omnibenchmark.git.cache.parse_repo_url", return_value="github.com/org/repo"
        ),
    ):
        # Should not raise — failure is logged, not re-raised
        _populate_git_cache(
            _mock_benchmark(repo_url="https://github.com/org/repo.git"),
            quiet=False,
            cores=1,
        )


@pytest.mark.short
def test_populate_git_cache_quiet_uses_progress_display(tmp_path):
    with (
        patch("omnibenchmark.cli.run.get_git_cache_dir", return_value=tmp_path),
        patch("omnibenchmark.cli.run.is_local_path", return_value=False),
        patch("omnibenchmark.cli.run.get_or_update_cached_repo"),
        patch(
            "omnibenchmark.git.cache.parse_repo_url", return_value="github.com/org/repo"
        ),
        patch("omnibenchmark.cli.progress.ProgressDisplay") as mock_pd,
    ):
        mock_pd.return_value = MagicMock()
        _populate_git_cache(
            _mock_benchmark(repo_url="https://github.com/org/repo.git"),
            quiet=True,
            cores=1,
        )
        mock_pd.return_value.start_task.assert_called_once()
        mock_pd.return_value.finish.assert_called_once()


@pytest.mark.short
def test_populate_git_cache_commit_in_cache_skips_fetch(tmp_path):
    """When commit already exists in the local dulwich repo, skip fetch."""
    full_commit = "a" * 40  # 40-char hex commit
    repo_cache_subdir = tmp_path / "github.com" / "org" / "repo"
    repo_cache_subdir.mkdir(parents=True)

    mock_repo = MagicMock()
    # repo[commit_bytes] succeeds → skip_fetch = True

    with (
        patch("omnibenchmark.cli.run.get_git_cache_dir", return_value=tmp_path),
        patch("omnibenchmark.cli.run.is_local_path", return_value=False),
        patch("omnibenchmark.cli.run.get_or_update_cached_repo") as mock_fetch,
        patch(
            "omnibenchmark.git.cache.parse_repo_url", return_value="github.com/org/repo"
        ),
        patch("dulwich.porcelain.open_repo", return_value=mock_repo),
    ):
        _populate_git_cache(
            _mock_benchmark(
                repo_url="https://github.com/org/repo.git", commit=full_commit
            ),
            quiet=False,
            cores=1,
        )
        # Fetch should NOT have been called since commit is in cache
        mock_fetch.assert_not_called()


@pytest.mark.short
def test_populate_git_cache_future_exception(tmp_path):
    """When future.result() raises, the error is captured in failed list."""
    from concurrent.futures import Future

    failing_future = Future()
    failing_future.set_exception(RuntimeError("unexpected failure"))

    with (
        patch("omnibenchmark.cli.run.get_git_cache_dir", return_value=tmp_path),
        patch("omnibenchmark.cli.run.is_local_path", return_value=False),
        patch("omnibenchmark.cli.run.get_or_update_cached_repo"),
        patch(
            "omnibenchmark.git.cache.parse_repo_url", return_value="github.com/org/repo"
        ),
        patch(
            "concurrent.futures.ThreadPoolExecutor.submit", return_value=failing_future
        ),
    ):
        # Should not raise — failure is logged
        _populate_git_cache(
            _mock_benchmark(repo_url="https://github.com/org/repo.git"),
            quiet=False,
            cores=1,
        )


@pytest.mark.short
def test_run_module_filter_flag_logs_warning():
    with patch("omnibenchmark.cli.run._run_benchmark") as mock_rb:
        runner = CliRunner()
        runner.invoke(run, ["tests/data/mock_benchmark.yaml", "-m", "M1"])
        mock_rb.assert_called_once()
        assert mock_rb.call_args.kwargs.get("module_filter") == "M1"


@pytest.mark.short
def test_populate_git_cache_metric_collector_repo(tmp_path):
    """metric_collectors with a repo URL are added to the repos dict."""
    mock_collector = MagicMock()
    mock_collector.repository = MagicMock()
    mock_collector.repository.url = "https://github.com/org/metrics.git"
    mock_collector.repository.commit = "abc1234"

    mock_b = _mock_benchmark()
    mock_b.model.metric_collectors = [mock_collector]

    with (
        patch("omnibenchmark.cli.run.get_git_cache_dir", return_value=tmp_path),
        patch("omnibenchmark.cli.run.is_local_path", return_value=False),
        patch("omnibenchmark.cli.run.get_or_update_cached_repo") as mock_fetch,
        patch(
            "omnibenchmark.git.cache.parse_repo_url",
            return_value="github.com/org/metrics",
        ),
    ):
        _populate_git_cache(mock_b, quiet=False, cores=1)
        mock_fetch.assert_called_once()


@pytest.mark.short
def test_populate_git_cache_commit_lookup_keyerror(tmp_path):
    """When dulwich raises KeyError for commit lookup, fall through to fetch."""
    full_commit = "b" * 40
    repo_cache_subdir = tmp_path / "github.com" / "org" / "repo"
    repo_cache_subdir.mkdir(parents=True)

    mock_repo = MagicMock()
    mock_repo.__getitem__ = MagicMock(side_effect=KeyError("not found"))

    with (
        patch("omnibenchmark.cli.run.get_git_cache_dir", return_value=tmp_path),
        patch("omnibenchmark.cli.run.is_local_path", return_value=False),
        patch("omnibenchmark.cli.run.get_or_update_cached_repo") as mock_fetch,
        patch(
            "omnibenchmark.git.cache.parse_repo_url", return_value="github.com/org/repo"
        ),
        patch("dulwich.porcelain.open_repo", return_value=mock_repo),
    ):
        _populate_git_cache(
            _mock_benchmark(
                repo_url="https://github.com/org/repo.git", commit=full_commit
            ),
            quiet=False,
            cores=1,
        )
        # KeyError → skip_fetch=False → fetch is called
        mock_fetch.assert_called_once()


@pytest.mark.short
def test_run_benchmark_remote_storage_true_bool_appended(tmp_path):
    """A True boolean value in storage args should add a bare --flag."""
    (tmp_path / "Snakefile").write_text("rule all: input: []")
    with (
        patch("omnibenchmark.cli.run.BenchmarkExecution") as mock_be,
        patch("omnibenchmark.cli.run._populate_git_cache"),
        patch("omnibenchmark.cli.run._generate_explicit_snakefile"),
        patch("omnibenchmark.cli.run.write_run_manifest"),
        patch("omnibenchmark.cli.run._run_snakemake") as mock_snakemake,
        patch("omnibenchmark.remote.storage.get_storage_from_benchmark"),
        patch(
            "omnibenchmark.remote.storage.remote_storage_snakemake_args"
        ) as mock_storage_args,
    ):
        mock_b = MagicMock()
        mock_b.get_benchmark_software_backend.return_value = SoftwareBackendEnum.host
        mock_be.return_value = mock_b
        mock_storage_args.return_value = {"use-conda": True}
        _run_benchmark(
            benchmark_path="tests/data/mock_benchmark.yaml",
            cores=1,
            dry=False,
            continue_on_error=False,
            out_dir=str(tmp_path),
            debug=False,
            dirty=False,
            use_remote_storage=True,
        )
        _, kwargs = mock_snakemake.call_args
        extra = kwargs["extra_snakemake_args"]
        assert "--use-conda" in extra
