import hashlib
import json
import os
import re
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

import pytest
import yaml

from omni.benchmark import Benchmark
from omni.io.MinIOStorage import MinIOStorage
from tests.cli.cli_setup import OmniCLISetup
from tests.io.MinIOStorage_setup import MinIOSetup, TmpMinIOStorage

if not sys.platform == "linux":
    pytest.skip(
        "for GHA, only works on linux (https://docs.github.com/en/actions/using-containerized-services/about-service-containers#about-service-containers)",
        allow_module_level=True,
    )
else:
    tempdir = Path(tempfile.gettempdir()) / "ob_test_benchmark004"
    # benchmark_data_path = Path("tests/data")
    # benchmark_path = benchmark_data_path / "Benchmark_004.yaml"
    benchmark_data = Path("..") / "data"
    benchmark_data_path = Path(__file__).parent / benchmark_data
    benchmark_path = str(benchmark_data_path / "Benchmark_004.yaml")

    with open(benchmark_path, "r") as fh:
        yaml.safe_load(fh)
        benchmark_obj = Benchmark(Path(benchmark_path))

    minio_testcontainer = MinIOSetup(sys.platform == "linux")


class TestCLIMinIOStorage:
    def test_archive_config(self, monkeypatch: pytest.MonkeyPatch):
        expected_output = """
        Created archive:
        """
        if not sys.platform == "linux":
            pytest.skip(
                "for GHA, only works on linux (https://docs.github.com/en/actions/using-containerized-services/about-service-containers#about-service-containers)",
                allow_module_level=True,
            )

        with TmpMinIOStorage(minio_testcontainer) as tmp:
            tmp.setup(in_dir=benchmark_data_path, out_dir=tempdir)
            monkeypatch.chdir(tmp.out_dir)
            subprocess.run(
                ["ob", "run", "benchmark", "--benchmark", tmp.benchmark_file],
                cwd=tmp.out_dir,
            )
            subprocess.run(
                ["ob", "storage", "create-version", "--benchmark", tmp.benchmark_file],
                cwd=tmp.out_dir,
            )
            with OmniCLISetup() as omni:
                result = omni.call(
                    [
                        "info",
                        "archive",
                        "--benchmark",
                        str(tmp.benchmark_file),
                    ]
                )
                assert result.exit_code == 0
                assert clean(result.output).startswith(clean(expected_output))
                outfile = f"{benchmark_obj.get_benchmark_name()}_{benchmark_obj.get_converter().get_version()}.zip"
                assert Path(outfile).exists()
                with zipfile.ZipFile(outfile, "r") as f:
                    files = f.namelist()
                assert "/".join(tmp.benchmark_file.parts[1:]) in files

    def test_archive_code(self, monkeypatch: pytest.MonkeyPatch):
        expected_output = """
        Created archive:
        """
        if not sys.platform == "linux":
            pytest.skip(
                "for GHA, only works on linux (https://docs.github.com/en/actions/using-containerized-services/about-service-containers#about-service-containers)",
                allow_module_level=True,
            )

        with TmpMinIOStorage(minio_testcontainer) as tmp:
            tmp.setup(in_dir=benchmark_data_path, out_dir=tempdir)
            monkeypatch.chdir(tmp.out_dir)
            subprocess.run(
                ["ob", "run", "benchmark", "--benchmark", tmp.benchmark_file],
                cwd=tmp.out_dir,
            )
            subprocess.run(
                ["ob", "storage", "create-version", "--benchmark", tmp.benchmark_file],
                cwd=tmp.out_dir,
            )
            with OmniCLISetup() as omni:
                result = omni.call(
                    ["info", "archive", "--benchmark", str(tmp.benchmark_file), "-c"]
                )
                assert result.exit_code == 0
                assert clean(result.output).startswith(clean(expected_output))
                outfile = f"{benchmark_obj.get_benchmark_name()}_{benchmark_obj.get_converter().get_version()}.zip"
                assert Path(outfile).exists()
                with zipfile.ZipFile(outfile, "r") as f:
                    files = f.namelist()
                print(str(tmp.benchmark_file))
                print(files)
                testfile = (
                    ".snakemake/repos/6f0b9a981e9a3c8329af1e67ce49631e/config.cfg"
                )
                assert testfile in files

    def test_archive_results(self, monkeypatch: pytest.MonkeyPatch):
        expected_output = """
        Created archive:
        """
        if not sys.platform == "linux":
            pytest.skip(
                "for GHA, only works on linux (https://docs.github.com/en/actions/using-containerized-services/about-service-containers#about-service-containers)",
                allow_module_level=True,
            )

        with TmpMinIOStorage(minio_testcontainer) as tmp:
            tmp.setup(in_dir=benchmark_data_path, out_dir=tempdir)
            monkeypatch.chdir(tmp.out_dir)
            subprocess.run(
                ["ob", "run", "benchmark", "--benchmark", tmp.benchmark_file],
                cwd=tmp.out_dir,
            )
            subprocess.run(
                ["ob", "storage", "create-version", "--benchmark", tmp.benchmark_file],
                cwd=tmp.out_dir,
            )
            with OmniCLISetup() as omni:
                result = omni.call(
                    ["info", "archive", "--benchmark", str(tmp.benchmark_file), "-r"]
                )
                assert result.exit_code == 0
                assert clean(result.output).startswith(clean(expected_output))
                outfile = f"{benchmark_obj.get_benchmark_name()}_{benchmark_obj.get_converter().get_version()}.zip"
                assert Path(outfile).exists()
                with zipfile.ZipFile(outfile, "r") as f:
                    files = f.namelist()
                print(str(tmp.benchmark_file))
                print(files)
                testfile = "out/data/D2/default/D2.meta.json"
                assert testfile in files


def clean(output: str) -> str:
    output = output.strip()
    output = output.replace("    ", "")

    # Replace different newline characters with a single '\n'
    normalized_output = re.sub(r"\r\n|\r", "\n", output)

    # Replace multiple spaces and tabs with a single space
    normalized_output = re.sub(r"[ \t]+", " ", normalized_output)

    return normalized_output


def cleanup_buckets_on_exit():
    """Cleanup a testing directory once we are finished."""
    TmpMinIOStorage(minio_testcontainer).cleanup_buckets()


@pytest.fixture(scope="session", autouse=True)
def cleanup(request):
    """Cleanup a testing directory once we are finished."""
    request.addfinalizer(cleanup_buckets_on_exit)
