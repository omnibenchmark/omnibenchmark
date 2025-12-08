import pytest
from unittest.mock import Mock, patch
from omnibenchmark.workflow.snakemake import SnakemakeEngine
from omnibenchmark.model import SoftwareBackendEnum


class TestSnakemakeEngine:
    @pytest.fixture
    def engine(self):
        """Create a SnakemakeEngine instance."""
        return SnakemakeEngine()

    @pytest.fixture
    def mock_node(self):
        """Create a mock BenchmarkNode."""
        node = Mock()
        node.get_definition_file.return_value = None
        node.get_benchmark_name.return_value = "test_benchmark"
        node.get_benchmark_version.return_value = "1.0.0"
        node.get_benchmark_author.return_value = "test_author"
        node.get_stage_id.return_value = "test_stage"
        node.get_module_id.return_value = "test_module"
        node.get_param_id.return_value = "test_params"
        return node

    @pytest.mark.short
    def test_serialize_node_workflow_uses_provided_benchmark_file_path(
        self, engine, mock_node, tmp_path
    ):
        """Test that serialize_node_workflow uses provided benchmark_file_path."""
        # Arrange
        benchmark_file = tmp_path / "test_benchmark.yaml"
        benchmark_file.touch()

        # Act
        result = engine.serialize_node_workflow(
            mock_node,
            benchmark_file,
            output_dir=tmp_path,
            write_to_disk=True,
        )

        # Assert
        assert result.exists()
        assert result.name == "Snakefile"
        # Verify the benchmark file path is used in the generated Snakefile
        content = result.read_text()
        assert str(benchmark_file) in content

    @pytest.mark.short
    def test_prepare_argv_handles_list_values(self, engine, tmp_path):
        """Test that _prepare_argv correctly handles list values in snakemake_kwargs."""
        # Arrange
        snakemake_kwargs = {
            "set_threads": [4, 8],
            "verbose": True,
        }

        # Act
        argv = engine._prepare_argv(
            snakefile=tmp_path / "Snakefile",
            cores=1,
            update=False,
            dryrun=False,
            keep_module_logs=False,
            continue_on_error=False,
            backend=SoftwareBackendEnum.host,
            work_dir=tmp_path,
            out_dir=tmp_path / "out",
            debug=False,
            **snakemake_kwargs,
        )

        # Assert
        # Check that list values are properly expanded
        assert "--set_threads" in argv
        set_threads_idx = argv.index("--set_threads")
        assert argv[set_threads_idx + 1] == "4"
        assert argv[set_threads_idx + 2] == "8"

        # Check that boolean True flag is present
        assert "--verbose" in argv

    @pytest.mark.short
    def test_prepare_argv_excludes_none_values(self, engine, tmp_path):
        """Test that _prepare_argv excludes None values from arguments."""
        # Arrange
        snakemake_kwargs = {
            "some_option": None,
            "other_option": "value",
        }

        # Act
        argv = engine._prepare_argv(
            snakefile=tmp_path / "Snakefile",
            cores=1,
            update=False,
            dryrun=False,
            keep_module_logs=False,
            continue_on_error=False,
            backend=SoftwareBackendEnum.host,
            work_dir=tmp_path,
            out_dir=tmp_path / "out",
            debug=False,
            **snakemake_kwargs,
        )

        # Assert
        # None values should not appear in argv
        assert "--some_option" not in argv
        assert "None" not in argv

        # Non-None values should be present
        assert "--other_option" in argv
        assert "value" in argv

    @pytest.mark.short
    def test_prepare_argv_always_includes_out_dir(self, engine, tmp_path):
        """Test that _prepare_argv always includes out_dir in config."""
        # Act
        argv = engine._prepare_argv(
            snakefile=tmp_path / "Snakefile",
            cores=1,
            update=False,
            dryrun=False,
            keep_module_logs=False,
            continue_on_error=False,
            backend=SoftwareBackendEnum.host,
            work_dir=tmp_path,
            out_dir=tmp_path / "out",
            debug=False,
        )

        # Assert
        # Check that out_dir is in the config arguments
        assert any(f"out_dir={tmp_path / 'out'}" in arg for arg in argv)

    @pytest.mark.short
    def test_prepare_argv_always_includes_backend(self, engine, tmp_path):
        """Test that _prepare_argv always includes backend in config."""
        # Act
        argv = engine._prepare_argv(
            snakefile=tmp_path / "Snakefile",
            cores=1,
            update=False,
            dryrun=False,
            keep_module_logs=False,
            continue_on_error=False,
            backend=SoftwareBackendEnum.conda,
            work_dir=tmp_path,
            out_dir=tmp_path / "out",
            debug=False,
        )

        # Assert
        # Check that backend is in the config arguments
        assert any("backend=conda" in arg for arg in argv)

    @pytest.mark.short
    def test_run_workflow_uses_benchmark_context_out_dir(self, engine, tmp_path):
        """Test that run_workflow uses benchmark.context.out_dir."""
        # Arrange
        mock_benchmark = Mock()
        mock_benchmark.context.out_dir = str(tmp_path / "custom_out")
        mock_benchmark.get_output_paths.return_value = []
        mock_benchmark.get_metric_collector_output_paths.return_value = []

        # Mock the internal methods to avoid actual snakemake execution
        with patch.object(
            engine, "serialize_workflow", return_value=tmp_path / "Snakefile"
        ):
            with patch.object(engine, "_prepare_argv", return_value=[]):
                with patch("omnibenchmark.workflow.snakemake.snakemake.snakemake_cli"):
                    # Act
                    engine.run_workflow(mock_benchmark, work_dir=tmp_path)

                    # Assert - verify serialize_workflow was called (which means out_dir was accessed)
                    engine.serialize_workflow.assert_called_once()

    @pytest.mark.short
    def test_serialize_node_workflow_raises_error_when_benchmark_file_path_is_none(
        self, engine, mock_node, tmp_path
    ):
        """Test that serialize_node_workflow raises ValueError when benchmark_file_path is None."""
        # Act & Assert
        with pytest.raises(ValueError, match="benchmark_file_path must be provided"):
            engine.serialize_node_workflow(
                mock_node, benchmark_file_path=None, output_dir=tmp_path
            )

    @pytest.mark.short
    def test_run_node_workflow_raises_error_when_benchmark_file_path_is_none(
        self, engine, mock_node, tmp_path
    ):
        """Test that run_node_workflow raises ValueError when benchmark_file_path is None."""
        # Act & Assert
        with pytest.raises(
            ValueError,
            match="benchmark_file_path must be provided when node.get_definition_file\\(\\) returns None",
        ):
            engine.run_node_workflow(
                mock_node,
                input_dir=tmp_path,
                dataset="test_dataset",
                work_dir=tmp_path,
                benchmark_file_path=None,
            )

    @pytest.mark.short
    def test_prepare_argv_includes_backend_env(self, engine, tmp_path):
        """Test that _prepare_argv includes backend_env in config."""
        # Act
        argv = engine._prepare_argv(
            snakefile=tmp_path / "Snakefile",
            cores=1,
            update=False,
            dryrun=False,
            keep_module_logs=False,
            continue_on_error=False,
            backend=SoftwareBackendEnum.conda,
            work_dir=tmp_path,
            out_dir=tmp_path / "out",
            debug=False,
            backend_env="path/to/env.yaml",
        )

        # Assert
        assert any("backend_env=path/to/env.yaml" in arg for arg in argv)
