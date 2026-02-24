"""
Test to verify run (with --module filter) function works correctly.
"""

import pytest
from pathlib import Path
from unittest.mock import patch

from click.testing import CliRunner
from omnibenchmark.cli.run import run


@pytest.mark.short
def test_run_module_basic_functionality(tmp_path):
    """Test that run with --module flag and --dry generates a snakefile."""

    test_data_dir = Path(__file__).parent.parent / "data"
    benchmark_file = test_data_dir / "mock_benchmark.yaml"

    with patch("omnibenchmark.cli.run._populate_git_cache"):
        with patch("omnibenchmark.cli.run._generate_explicit_snakefile") as mock_gen:
            mock_gen.return_value = []
            with patch("omnibenchmark.cli.run.write_run_manifest"):
                runner = CliRunner()
                runner.invoke(
                    run,
                    [str(benchmark_file), "--module", "D1", "--dry"],
                    catch_exceptions=False,
                )

                # Verify BenchmarkExecution was called
                mock_gen.assert_called_once()
                args, kwargs = mock_gen.call_args
                # module_filter should be passed through
                assert kwargs.get("module_filter") == "D1" or "D1" in str(args)
