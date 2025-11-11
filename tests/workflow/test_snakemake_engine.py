"""Tests for SnakemakeEngine coverage gaps."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from omnibenchmark.workflow.snakemake.snakemake import (
    SnakemakeEngine,
    _print_ob_flags_help,
)


class TestPrintObFlagsHelp:
    """Tests for the _print_ob_flags_help function."""

    @pytest.mark.short
    def test_print_ob_flags_help_benchmark(self, capsys):
        _print_ob_flags_help("benchmark")
        captured = capsys.readouterr()
        assert "HINT:" in captured.err
        assert "ob run benchmark --help" in captured.err

    @pytest.mark.short
    def test_print_ob_flags_help_module(self, capsys):
        _print_ob_flags_help("module")
        captured = capsys.readouterr()
        assert "ob run module --help" in captured.err


class TestPrepareArgv:
    """Tests for the _prepare_argv method."""

    @pytest.mark.short
    def test_prepare_argv_with_list_values(self):
        from omnibenchmark.model import SoftwareBackendEnum

        argv = SnakemakeEngine._prepare_argv(
            snakefile=Path("/tmp/Snakefile"),
            cores=1,
            update=False,
            dryrun=False,
            keep_module_logs=False,
            continue_on_error=False,
            backend=SoftwareBackendEnum.host,
            work_dir=Path("/tmp"),
            out_dir="/tmp/out",
            config=["key1=val1", "key2=val2"],
        )
        assert "--config" in argv
        assert "key1=val1" in argv
        assert "key2=val2" in argv

    @pytest.mark.short
    def test_prepare_argv_with_none_values(self):
        from omnibenchmark.model import SoftwareBackendEnum

        argv = SnakemakeEngine._prepare_argv(
            snakefile=Path("/tmp/Snakefile"),
            cores=1,
            update=False,
            dryrun=False,
            keep_module_logs=False,
            continue_on_error=False,
            backend=SoftwareBackendEnum.host,
            work_dir=Path("/tmp"),
            out_dir="/tmp/out",
            some_none_param=None,
        )
        assert "--some-none-param" not in argv


class TestSystemExitHandling:
    """Tests for SystemExit handling in run methods."""

    @pytest.mark.short
    @patch("omnibenchmark.workflow.snakemake.snakemake.parse_args")
    def test_run_workflow_system_exit(self, mock_parse_args, capsys):
        mock_parse_args.side_effect = SystemExit(2)

        engine = SnakemakeEngine()
        mock_benchmark = MagicMock()
        mock_benchmark.context.out_dir = Path("/tmp/out")
        mock_benchmark.get_nodes.return_value = []

        with patch.object(
            engine, "serialize_workflow", return_value=Path("/tmp/Snakefile")
        ):
            with pytest.raises(SystemExit):
                engine.run_workflow(mock_benchmark)

        captured = capsys.readouterr()
        assert "HINT:" in captured.err

    @pytest.mark.short
    @patch("omnibenchmark.workflow.snakemake.snakemake.parse_args")
    def test_run_node_workflow_system_exit(self, mock_parse_args, capsys, tmp_path):
        mock_parse_args.side_effect = SystemExit(2)

        engine = SnakemakeEngine()
        mock_node = MagicMock()
        mock_node.get_definition_file.return_value = Path("benchmark.yaml")
        mock_node.get_benchmark_name.return_value = "test"
        mock_node.get_benchmark_version.return_value = "1.0"
        mock_node.get_benchmark_author.return_value = "author"

        with patch.object(
            engine, "serialize_node_workflow", return_value=Path("/tmp/Snakefile")
        ):
            with pytest.raises(SystemExit):
                engine.run_node_workflow(
                    mock_node,
                    input_dir=tmp_path,
                    dataset="test_dataset",
                    work_dir=tmp_path,
                )

        captured = capsys.readouterr()
        assert "HINT:" in captured.err

    @pytest.mark.short
    @patch("omnibenchmark.workflow.snakemake.snakemake.parse_args")
    def test_run_workflow_system_exit_zero_no_hint(self, mock_parse_args, capsys):
        """SystemExit(0) should not print the hint."""
        mock_parse_args.side_effect = SystemExit(0)

        engine = SnakemakeEngine()
        mock_benchmark = MagicMock()
        mock_benchmark.context.out_dir = Path("/tmp/out")

        with patch.object(
            engine, "serialize_workflow", return_value=Path("/tmp/Snakefile")
        ):
            with pytest.raises(SystemExit):
                engine.run_workflow(mock_benchmark)

        captured = capsys.readouterr()
        assert "HINT:" not in captured.err


class TestSerializeNodeWorkflow:
    """Tests for serialize_node_workflow edge cases."""

    @pytest.mark.short
    def test_serialize_node_workflow_no_benchmark_file(self, tmp_path):
        engine = SnakemakeEngine()
        mock_node = MagicMock()
        mock_node.get_definition_file.return_value = None

        with pytest.raises(ValueError, match="benchmark_file_path must be provided"):
            engine.serialize_node_workflow(
                mock_node,
                output_dir=tmp_path,
                benchmark_file_path=None,
            )
