"""Tests for parse_performance.py script"""

import tempfile
import csv
from pathlib import Path
import importlib.util

import pytest


# Import the module directly to avoid snakemake dependency issues
def _import_parse_performance():
    """Import parse_performance.py directly without going through package init"""
    script_path = (
        Path(__file__).parent.parent.parent
        / "omnibenchmark"
        / "workflow"
        / "snakemake"
        / "scripts"
        / "parse_performance.py"
    )
    spec = importlib.util.spec_from_file_location("parse_performance", script_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


parse_performance = _import_parse_performance()
write_combined_performance_file = parse_performance.write_combined_performance_file
combine_performances = parse_performance.combine_performances
read_performance = parse_performance.read_performance


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files"""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def sample_performance_file(temp_dir):
    """Create a sample performance file"""
    perf_file = temp_dir / "test_performance.txt"
    with open(perf_file, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["s", "max_rss", "max_vms", "io_in", "io_out"], delimiter="\t"
        )
        writer.writeheader()
        writer.writerow(
            {
                "s": "10.5",
                "max_rss": "1024",
                "max_vms": "2048",
                "io_in": "100",
                "io_out": "200",
            }
        )
    return perf_file


def test_read_performance(sample_performance_file):
    """Test reading a performance file"""
    records = list(read_performance(str(sample_performance_file)))

    assert len(records) == 1
    record = records[0]

    assert record["s"] == 10.5
    assert record["max_rss"] == 1024.0
    assert record["max_vms"] == 2048.0
    assert record["io_in"] == 100.0
    assert record["io_out"] == 200.0


def test_read_performance_with_na_values(temp_dir):
    """Test reading a performance file with NA values"""
    perf_file = temp_dir / "test_na_performance.txt"
    with open(perf_file, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["s", "max_rss", "max_vms"], delimiter="\t"
        )
        writer.writeheader()
        writer.writerow({"s": "NA", "max_rss": "1024", "max_vms": "none"})

    records = list(read_performance(str(perf_file)))

    assert len(records) == 1
    record = records[0]

    assert record["s"] == 0  # NA converted to 0
    assert record["max_rss"] == 1024.0
    assert record["max_vms"] == 0  # 'none' converted to 0


def test_combine_performances_empty_list(temp_dir):
    """Test combining with empty file list"""
    result = combine_performances(temp_dir, [])
    assert result == []


def test_combine_performances_single_file(temp_dir, sample_performance_file):
    """Test combining a single performance file"""
    # Create directory structure for module extraction
    module_dir = (
        temp_dir / "data" / "dataset" / "default" / "method" / "module1" / "params"
    )
    module_dir.mkdir(parents=True)
    perf_in_module = module_dir / "test_performance.txt"

    # Copy sample file to module directory
    with open(sample_performance_file, "r") as src:
        with open(perf_in_module, "w") as dst:
            dst.write(src.read())

    result = combine_performances(temp_dir, [str(perf_in_module)])

    assert len(result) == 1
    record = result[0]

    assert record["s"] == 10.5
    assert record["module"] == "module1"
    assert record["path"] == str(perf_in_module)
    assert "params" in record


def test_combine_performances_multiple_files(temp_dir):
    """Test combining multiple performance files"""
    # Create two module directories with performance files
    for i in range(2):
        module_dir = (
            temp_dir
            / "data"
            / "dataset"
            / "default"
            / "method"
            / f"module{i}"
            / "params"
        )
        module_dir.mkdir(parents=True)
        perf_file = module_dir / "test_performance.txt"

        with open(perf_file, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["s", "max_rss"], delimiter="\t")
            writer.writeheader()
            writer.writerow({"s": str(10 + i), "max_rss": str(1000 * (i + 1))})

    # Get all performance files
    perf_files = list(temp_dir.glob("**/test_performance.txt"))
    perf_files_str = [str(f) for f in perf_files]

    result = combine_performances(temp_dir, perf_files_str)

    assert len(result) == 2
    assert result[0]["s"] in [10.0, 11.0]
    assert result[1]["s"] in [10.0, 11.0]


def test_write_combined_performance_file_empty(temp_dir):
    """Test writing combined performance file with no input files"""
    write_combined_performance_file(temp_dir, [])

    output_file = temp_dir / "performances.tsv"
    # An empty file should be created to satisfy Snakemake output requirements
    assert output_file.exists()
    # The file should be empty (no headers, no content)
    assert output_file.stat().st_size == 0


def test_write_combined_performance_file(temp_dir):
    """Test writing combined performance file"""
    # Create a performance file
    module_dir = (
        temp_dir / "data" / "dataset" / "default" / "method" / "testmodule" / "params"
    )
    module_dir.mkdir(parents=True)
    perf_file = module_dir / "test_performance.txt"

    with open(perf_file, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["s", "max_rss", "max_vms"], delimiter="\t"
        )
        writer.writeheader()
        writer.writerow({"s": "15.5", "max_rss": "2048", "max_vms": "4096"})

    write_combined_performance_file(temp_dir, [str(perf_file)])

    output_file = temp_dir / "performances.tsv"
    assert output_file.exists()

    # Read and verify the output
    with open(output_file, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)

    assert len(rows) == 1
    row = rows[0]

    assert row["s"] == "15.5"
    assert row["max_rss"] == "2048.0"
    assert row["max_vms"] == "4096.0"
    assert row["module"] == "testmodule"
    assert row["path"] == str(perf_file)


def test_write_combined_performance_file_creates_directory(temp_dir):
    """Test that write_combined_performance_file creates output directory if it doesn't exist"""
    output_dir = temp_dir / "new_output_dir"
    assert not output_dir.exists()

    write_combined_performance_file(output_dir, [])

    assert output_dir.exists()
