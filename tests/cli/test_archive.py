import zipfile
from pathlib import Path
from typing import Optional


from omnibenchmark.benchmark import Benchmark
from tests.cli.cli_setup import OmniCLISetup

from .fixtures import minio_storage, _minio_container  # noqa: F401


def do_first_run(clisetup, file: str, cwd: Optional[str] = None):
    run1 = clisetup.call(
        [
            "run",
            "benchmark",
            "--benchmark",
            file,
        ],
        cwd=cwd,
    )
    assert run1.returncode == 0


def test_archive_config(minio_storage):  # noqa: F811
    """Test archiving a benchmark configuration."""
    expected_output = "Created archive:"

    with OmniCLISetup() as omni:
        # First run the benchmark to generate data
        do_first_run(omni, str(minio_storage.benchmark_file), cwd=minio_storage.out_dir)

        # Create a version which is needed before archiving
        run2 = omni.call(
            [
                "storage",
                "create-version",
                "--benchmark",
                str(minio_storage.benchmark_file),
            ],
            cwd=minio_storage.out_dir,
        )
        assert run2.returncode == 0

        # Now run the archive command
        run3 = omni.call(
            [
                "storage",
                "archive",
                "--benchmark",
                str(minio_storage.benchmark_file),
            ],
            cwd=minio_storage.out_dir,
        )

        # Check the results
        assert run3.returncode == 0
        assert run3.stdout.startswith(expected_output)

        # Get the benchmark name and version to check the archive file
        benchmark = Benchmark(Path(minio_storage.benchmark_file))
        outfile = f"{benchmark.get_benchmark_name()}_{benchmark.get_converter().get_version()}.zip"

        # Check that the archive exists and contains the benchmark file
        archive_path = Path(minio_storage.out_dir) / outfile
        assert archive_path.exists()

        with zipfile.ZipFile(archive_path, "r") as f:
            files = f.namelist()

        # Check that the benchmark file is in the archive
        relative_benchmark_path = "/".join(Path(minio_storage.benchmark_file).parts[1:])
        assert relative_benchmark_path in files


def test_archive_code(minio_storage):  # noqa: F811
    """Test archiving a benchmark with code files."""
    expected_output = "Created archive:"

    with OmniCLISetup() as omni:
        # First run the benchmark to generate data
        do_first_run(omni, str(minio_storage.benchmark_file), cwd=minio_storage.out_dir)

        # Create a version which is needed before archiving
        run2 = omni.call(
            [
                "storage",
                "create-version",
                "--benchmark",
                str(minio_storage.benchmark_file),
            ],
            cwd=minio_storage.out_dir,
        )
        assert run2.returncode == 0

        # Now run the archive command with -c to include code
        run3 = omni.call(
            [
                "storage",
                "archive",
                "--benchmark",
                str(minio_storage.benchmark_file),
                "-c",  # Include code
            ],
            cwd=minio_storage.out_dir,
        )

        # Check the results
        assert run3.returncode == 0
        assert run3.stdout.startswith(expected_output)

        # Get the benchmark name and version to check the archive file
        benchmark = Benchmark(Path(minio_storage.benchmark_file))
        outfile = f"{benchmark.get_benchmark_name()}_{benchmark.get_converter().get_version()}.zip"

        # Check that the archive exists and contains the benchmark file
        archive_path = Path(minio_storage.out_dir) / outfile
        assert archive_path.exists()

        with zipfile.ZipFile(archive_path, "r") as f:
            files = f.namelist()

        # Check for a specific code file in the archive
        testfile = ".snakemake/repos/cc6e274b72f3ef949eae61aaf7ee3b7653f0bbf92aa2ea879a2359acb137dbb4/config.cfg"
        assert testfile in files


def test_archive_results(minio_storage):  # noqa: F811
    """Test archiving benchmark results."""
    expected_output = "Created archive:"

    with OmniCLISetup() as omni:
        # First run the benchmark to generate data
        do_first_run(omni, str(minio_storage.benchmark_file), cwd=minio_storage.out_dir)

        # Create a version which is needed before archiving
        run2 = omni.call(
            [
                "storage",
                "create-version",
                "--benchmark",
                str(minio_storage.benchmark_file),
            ],
            cwd=minio_storage.out_dir,
        )
        assert run2.returncode == 0

        # Now run the archive command with -r to include results
        run3 = omni.call(
            [
                "storage",
                "archive",
                "--benchmark",
                str(minio_storage.benchmark_file),
                "-r",  # Include results
            ],
            cwd=minio_storage.out_dir,
        )

        # Check the results
        #
        assert run3.returncode == 0
        assert run3.stdout.startswith(expected_output)

        # Get the benchmark name and version to check the archive file
        benchmark = Benchmark(Path(minio_storage.benchmark_file))
        outfile = f"{benchmark.get_benchmark_name()}_{benchmark.get_converter().get_version()}.zip"

        # Check that the archive exists and contains the benchmark file
        archive_path = Path(minio_storage.out_dir) / outfile
        assert archive_path.exists()

        with zipfile.ZipFile(archive_path, "r") as f:
            files = f.namelist()

        # Check for a specific results file in the archive
        testfile = "out/data/D2/default/D2.meta.json"
        assert testfile in files


def test_archive_compression(minio_storage):  # noqa: F811
    """Test archiving with different compression."""
    expected_output = "Created archive:"

    with OmniCLISetup() as omni:
        # First run the benchmark to generate data
        do_first_run(omni, str(minio_storage.benchmark_file), cwd=minio_storage.out_dir)

        # Create a version which is needed before archiving
        run2 = omni.call(
            [
                "storage",
                "create-version",
                "--benchmark",
                str(minio_storage.benchmark_file),
            ],
            cwd=minio_storage.out_dir,
        )
        assert run2.returncode == 0

        # Now run the archive command with -r and --compression
        run3 = omni.call(
            [
                "storage",
                "archive",
                "--benchmark",
                str(minio_storage.benchmark_file),
                "-r",  # Include results
                "--compression",
                "bzip2",
            ],
            cwd=minio_storage.out_dir,
        )

        # Check the results
        assert run3.returncode == 0
        assert run3.stdout.startswith(expected_output)

        # Get the benchmark name and version to check the archive file
        benchmark = Benchmark(Path(minio_storage.benchmark_file))
        outfile = f"{benchmark.get_benchmark_name()}_{benchmark.get_converter().get_version()}.bz2"

        # Check that the archive exists and contains the benchmark file
        archive_path = Path(minio_storage.out_dir) / outfile
        assert archive_path.exists()

        with zipfile.ZipFile(archive_path, "r") as f:
            files = f.namelist()

        # Check for a specific results file in the archive
        testfile = "out/data/D2/default/D2.meta.json"
        assert testfile in files
