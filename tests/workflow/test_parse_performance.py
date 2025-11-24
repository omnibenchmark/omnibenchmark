"""Unit tests for parse_performance.py script.

Tests the performance file parsing and combining functionality.
"""

import pytest
from pathlib import Path
import csv

from omnibenchmark.workflow.snakemake.scripts.parse_performance import (
    read_performance,
    tokenize,
    read_params,
    combine_performances,
    write_combined_performance_file,
)


@pytest.mark.short
class TestReadPerformance:
    """Test reading individual performance files."""

    def test_read_performance_valid_file(self, tmp_path):
        """Test reading a valid performance file."""
        perf_file = tmp_path / "test_performance.txt"
        perf_file.write_text(
            "s\tmax_rss\tmax_vms\tmax_uss\tmax_pss\tio_in\tio_out\tmean_load\tcpu_time\th:m:s\n"
            "10.5\t1024\t2048\t512\t768\t100\t200\t0.5\t8.2\t00:00:10\n"
        )

        records = list(read_performance(str(perf_file)))

        assert len(records) == 1
        record = records[0]

        # Check that h:m:s is removed
        assert "h:m:s" not in record

        # Check values are converted to float
        assert record["s"] == 10.5
        assert record["max_rss"] == 1024.0
        assert record["cpu_time"] == 8.2

    def test_read_performance_handles_na_values(self, tmp_path):
        """Test that NA values are converted to 0."""
        perf_file = tmp_path / "test_performance.txt"
        perf_file.write_text("s\tmax_rss\tmax_vms\n" "10.5\tNA\tnone\n" "5.2\t\tnull\n")

        records = list(read_performance(str(perf_file)))

        assert len(records) == 2

        # First record
        assert records[0]["s"] == 10.5
        assert records[0]["max_rss"] == 0
        assert records[0]["max_vms"] == 0

        # Second record
        assert records[1]["s"] == 5.2
        assert records[1]["max_rss"] == 0
        assert records[1]["max_vms"] == 0

    def test_read_performance_handles_invalid_values(self, tmp_path):
        """Test that invalid numeric values are converted to 0."""
        perf_file = tmp_path / "test_performance.txt"
        perf_file.write_text("s\tmax_rss\tcomment\n" "10.5\t1024\tsome_text\n")

        records = list(read_performance(str(perf_file)))

        assert len(records) == 1
        assert records[0]["s"] == 10.5
        assert records[0]["max_rss"] == 1024.0
        assert records[0]["comment"] == 0  # Invalid text converted to 0


@pytest.mark.short
class TestTokenize:
    """Test tokenizing file paths into stage/module/params triples."""

    def test_tokenize_valid_path(self):
        """Test tokenizing a valid output path."""
        output_path = Path("out")
        file_path = "out/stage1/module1/params1/stage2/module2/params2/file.txt"

        triples = tokenize(output_path, file_path)

        assert len(triples) == 2
        assert triples[0] == ("stage1", "module1", "params1")
        assert triples[1] == ("stage2", "module2", "params2")

    def test_tokenize_single_triple(self):
        """Test tokenizing a path with single stage/module/params."""
        output_path = Path("output")
        file_path = "output/data/generator/.abc123/result.txt"

        triples = tokenize(output_path, file_path)

        assert len(triples) == 1
        assert triples[0] == ("data", "generator", ".abc123")

    def test_tokenize_invalid_path(self):
        """Test tokenizing a path that doesn't contain the output dir."""
        output_path = Path("out")
        file_path = "other/path/file.txt"

        triples = tokenize(output_path, file_path)

        assert triples == []


@pytest.mark.short
class TestReadParams:
    """Test reading parameter files."""

    def test_read_params_with_parameters(self, tmp_path):
        """Test reading parameters from a structured output directory."""
        output_path = tmp_path / "out"
        output_path.mkdir()

        # Create directory structure: out/stage1/module1/params1/
        stage_dir = output_path / "stage1" / "module1" / ".abc123"
        stage_dir.mkdir(parents=True)

        # Create parameters.json
        param_file = stage_dir / "parameters.json"
        param_file.write_text('{"option": "value1"}')

        # Test file path
        file_path = str(output_path / "stage1" / "module1" / ".abc123" / "result.txt")

        params = read_params(output_path, file_path)

        assert "stage1" in params
        assert "module1" in params
        assert ".abc123" in params
        assert '{"option": "value1"}' in params

    def test_read_params_skips_default(self, tmp_path):
        """Test that parameters with 'default' in the name are skipped."""
        output_path = tmp_path / "out"
        output_path.mkdir()

        # Create directory structure with 'default' params
        stage_dir = output_path / "stage1" / "module1" / "default"
        stage_dir.mkdir(parents=True)

        param_file = stage_dir / "parameters.json"
        param_file.write_text('{"option": "value1"}')

        file_path = str(output_path / "stage1" / "module1" / "default" / "result.txt")

        params = read_params(output_path, file_path)

        # Should return empty string since 'default' is skipped
        assert params == ""


@pytest.mark.short
class TestCombinePerformances:
    """Test combining multiple performance files."""

    def test_combine_performances(self, tmp_path):
        """Test combining performance data from multiple files."""
        # Create output directory structure
        output_path = tmp_path / "out"
        output_path.mkdir()

        # Create first performance file
        stage1_dir = output_path / "stage1" / "module1" / ".abc123"
        stage1_dir.mkdir(parents=True)
        perf1 = stage1_dir / "perf_performance.txt"
        perf1.write_text("s\tmax_rss\n10.5\t1024\n")

        # Create parameters.json for first file
        (stage1_dir / "parameters.json").write_text('{"p1": "v1"}')

        # Create second performance file
        stage2_dir = output_path / "stage2" / "module2" / ".def456"
        stage2_dir.mkdir(parents=True)
        perf2 = stage2_dir / "perf_performance.txt"
        perf2.write_text("s\tmax_rss\n5.2\t512\n")

        # Create parameters.json for second file
        (stage2_dir / "parameters.json").write_text('{"p2": "v2"}')

        performance_files = [str(perf1), str(perf2)]

        rows = combine_performances(output_path, performance_files)

        assert len(rows) == 2

        # Check first row
        assert rows[0]["s"] == 10.5
        assert rows[0]["max_rss"] == 1024.0
        assert rows[0]["module"] == "module1"
        assert rows[0]["path"] == str(perf1)
        assert "stage1" in rows[0]["params"]

        # Check second row
        assert rows[1]["s"] == 5.2
        assert rows[1]["max_rss"] == 512.0
        assert rows[1]["module"] == "module2"
        assert rows[1]["path"] == str(perf2)
        assert "stage2" in rows[1]["params"]

    def test_combine_performances_skips_missing_files(self, tmp_path):
        """Test that missing performance files are skipped."""
        output_path = tmp_path / "out"
        output_path.mkdir()

        performance_files = [
            str(tmp_path / "nonexistent.txt"),
            str(tmp_path / "also_missing.txt"),
        ]

        rows = combine_performances(output_path, performance_files)

        assert rows == []


@pytest.mark.short
class TestWriteCombinedPerformanceFile:
    """Test writing combined performance file."""

    def test_write_combined_performance_file(self, tmp_path):
        """Test writing the combined performances.tsv file."""
        output_path = tmp_path / "out"
        output_path.mkdir()

        # Create a performance file
        stage_dir = output_path / "stage1" / "module1" / ".abc123"
        stage_dir.mkdir(parents=True)
        perf_file = stage_dir / "test_performance.txt"
        perf_file.write_text("s\tmax_rss\n10.5\t1024\n")

        # Create parameters.json
        (stage_dir / "parameters.json").write_text('{"p": "v"}')

        write_combined_performance_file(output_path, [str(perf_file)])

        # Check that performances.tsv was created
        output_file = output_path / "performances.tsv"
        assert output_file.exists()

        # Read and verify the TSV content
        with open(output_file, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)

            assert len(rows) == 1
            assert float(rows[0]["s"]) == 10.5
            assert float(rows[0]["max_rss"]) == 1024.0
            assert rows[0]["module"] == "module1"

    def test_write_creates_output_dir(self, tmp_path):
        """Test that write_combined_performance_file creates the output directory."""
        output_path = tmp_path / "out"
        # Don't create the directory

        # This should create the directory
        write_combined_performance_file(output_path, [])

        assert output_path.exists()
        assert output_path.is_dir()

    def test_write_handles_empty_performance_list(self, tmp_path):
        """Test that writing with no performance files doesn't create the TSV."""
        output_path = tmp_path / "out"
        output_path.mkdir()

        write_combined_performance_file(output_path, [])

        # performances.tsv should not be created if there are no rows
        output_file = output_path / "performances.tsv"
        assert not output_file.exists()
