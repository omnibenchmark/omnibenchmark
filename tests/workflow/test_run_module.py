"""
Test to verify run_module function works correctly.
"""

import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock

from click.testing import CliRunner
from omnibenchmark.cli.run import run


@pytest.mark.short
def test_run_module_basic_functionality():
    """Test that run_module function creates proper temp directory."""

    # Mock BenchmarkExecution to avoid actual file operations
    with patch("omnibenchmark.cli.run.BenchmarkExecution") as mock_benchmark_exec:
        mock_benchmark = MagicMock()
        mock_nodes = [MagicMock()]
        mock_nodes[0].is_entrypoint.return_value = True
        mock_benchmark.get_nodes_by_module_id.return_value = mock_nodes
        mock_benchmark.get_benchmark_datasets.return_value = ["D1"]
        mock_benchmark_exec.return_value = mock_benchmark

        # Use test data
        test_data_dir = Path(__file__).parent.parent / "data"
        benchmark_file = test_data_dir / "mock_benchmark.yaml"

        # Use CliRunner to test the command
        runner = CliRunner()
        runner.invoke(run, [str(benchmark_file), "--module", "D1", "--dry"])

        # Verify BenchmarkExecution was called with proper arguments
        mock_benchmark_exec.assert_called_once()
        args, kwargs = mock_benchmark_exec.call_args

        # Verify the temp directory argument is a valid Path object
        temp_dir_arg = args[1]
        assert isinstance(temp_dir_arg, Path)
        assert temp_dir_arg.is_absolute()
