import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import pytest
import yaml

from omni.benchmark import Benchmark
from omni.io.MinIOStorage import MinIOStorage
from tests.cli.cli_setup import OmniCLISetup
from tests.io.MinIOStorage_setup import MinIOSetup, TmpMinIOStorage, create_remote_test

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

    benchmark_path = create_remote_test(
        minio_testcontainer, in_dir=benchmark_data_path, out_dir=tempdir
    )


class TestCLIMinIOStorage:
    def test_create_version(self):
        expected_output = """
        Create a new benchmark version
        """
        if not sys.platform == "linux":
            pytest.skip(
                "for GHA, only works on linux (https://docs.github.com/en/actions/using-containerized-services/about-service-containers#about-service-containers)",
                allow_module_level=True,
            )

        with TmpMinIOStorage(minio_testcontainer) as tmp:
            time.sleep(2)
            # TODO: to setup bucket
            _ = MinIOStorage(auth_options=tmp.auth_options, benchmark="benchmark004")

            subprocess.run(["ob", "run", "benchmark", "--benchmark", benchmark_path])
            with OmniCLISetup() as omni:
                result = omni.call(
                    [
                        "storage",
                        "create-version",
                        "--benchmark",
                        str(benchmark_path),
                    ]
                )
                assert result.exit_code == 0
                assert clean(result.output).startswith(clean(expected_output))
                ss = MinIOStorage(
                    auth_options=tmp.auth_options, benchmark="benchmark004"
                )
                ss.set_version("1.0")
                ss._get_objects()
                assert "versions/1.0.csv" in ss.files.keys()

    def test_list_files(self):
        expected_output = (
            """1a92ffdf9b074c131af83cd288b65e5f out/data/D1/default/D1.meta.json"""
        )
        if not sys.platform == "linux":
            pytest.skip(
                "for GHA, only works on linux (https://docs.github.com/en/actions/using-containerized-services/about-service-containers#about-service-containers)",
                allow_module_level=True,
            )

        with TmpMinIOStorage(minio_testcontainer) as tmp:
            time.sleep(2)
            # TODO: to setup bucket
            _ = MinIOStorage(auth_options=tmp.auth_options, benchmark="benchmark004")

            subprocess.run(["ob", "run", "benchmark", "--benchmark", benchmark_path])
            subprocess.run(
                ["ob", "storage", "create-version", "--benchmark", benchmark_path]
            )
            with OmniCLISetup() as omni:
                result = omni.call(
                    [
                        "storage",
                        "list",
                        "--benchmark",
                        str(benchmark_path),
                    ]
                )
                assert result.exit_code == 0
            assert clean(result.output).startswith(clean(expected_output))

    def test_download_files(self):
        expected_output = """Downloading 8 files with a total size of 2.7KiB"""
        if not sys.platform == "linux":
            pytest.skip(
                "for GHA, only works on linux (https://docs.github.com/en/actions/using-containerized-services/about-service-containers#about-service-containers)",
                allow_module_level=True,
            )

        with TmpMinIOStorage(minio_testcontainer) as tmp:
            tmp = TmpMinIOStorage(minio_testcontainer)
            time.sleep(2)
            # TODO: to setup bucket
            ss = MinIOStorage(auth_options=tmp.auth_options, benchmark="benchmark004")

            subprocess.run(["ob", "run", "benchmark", "--benchmark", benchmark_path])
            subprocess.run(
                ["ob", "storage", "create-version", "--benchmark", benchmark_path]
            )
            with OmniCLISetup() as omni:
                omni = OmniCLISetup()
                with omni.runner.isolated_filesystem(tempdir):
                    result = omni.call(
                        [
                            "storage",
                            "download",
                            "--benchmark",
                            str(benchmark_path),
                        ]
                    )
                    curtempdir = os.getcwd()
                assert result.exit_code == 0
            assert clean(result.output).startswith(clean(expected_output))
            ss.set_version("1.0")
            ss._get_objects()
            files = list(ss.files.keys())
            files = [f for f in files if Path(f).parents[-2].name == "out"]
            assert all(
                [Path(curtempdir) / f in list(tempdir.rglob("*")) for f in files]
            )
            shutil.rmtree(curtempdir)

        with TmpMinIOStorage(minio_testcontainer) as tmp:
            tmp = TmpMinIOStorage(minio_testcontainer)
            time.sleep(2)
            # TODO: to setup bucket
            ss = MinIOStorage(auth_options=tmp.auth_options, benchmark="benchmark004")

    def test_missing_S3_storage_credentials_1(self):
        expected_output = """Invalid S3 config. Missing access_key and secret_key in environment variables (OB_STORAGE_S3_ACCESS_KEY, OB_STORAGE_S3_SECRET_KEY) or OB_STORAGE_S3_CONFIG"""
        if not sys.platform == "linux":
            pytest.skip(
                "for GHA, only works on linux (https://docs.github.com/en/actions/using-containerized-services/about-service-containers#about-service-containers)",
                allow_module_level=True,
            )
        with TmpMinIOStorage(minio_testcontainer) as tmp:
            # os.environ["OB_STORAGE_S3_ACCESS_KEY"] = ""
            os.environ.pop("OB_STORAGE_S3_ACCESS_KEY", None)
            os.environ["OB_STORAGE_S3_SECRET_KEY"] = tmp.auth_options["secret_key"]
            _ = MinIOStorage(auth_options=tmp.auth_options, benchmark="benchmark004")

            subprocess.run(["ob", "run", "benchmark", "--benchmark", benchmark_path])
            with OmniCLISetup() as omni:
                result = omni.call(
                    [
                        "storage",
                        "create-version",
                        "--benchmark",
                        str(benchmark_path),
                    ]
                )
                assert result.exit_code == 1
                assert clean(result.output).startswith(clean(expected_output))

    def test_S3_storage_credentials_from_file(self):
        expected_output = """Create a new benchmark version"""
        if not sys.platform == "linux":
            pytest.skip(
                "for GHA, only works on linux (https://docs.github.com/en/actions/using-containerized-services/about-service-containers#about-service-containers)",
                allow_module_level=True,
            )
        with TmpMinIOStorage(minio_testcontainer) as tmp:
            storage_s3_json = benchmark_data_path / "storage_s3.json"
            json.dump(tmp.auth_options, storage_s3_json.open("w"))

            os.environ["OB_STORAGE_S3_CONFIG"] = str(storage_s3_json)
            os.environ.pop("OB_STORAGE_S3_ACCESS_KEY", None)
            os.environ.pop("OB_STORAGE_S3_SECRET_KEY", None)

            _ = MinIOStorage(auth_options=tmp.auth_options, benchmark="benchmark004")

            subprocess.run(["ob", "run", "benchmark", "--benchmark", benchmark_path])
            with OmniCLISetup() as omni:
                result = omni.call(
                    [
                        "storage",
                        "create-version",
                        "--benchmark",
                        str(benchmark_path),
                    ]
                )
                assert result.exit_code == 0
                assert clean(result.output).startswith(clean(expected_output))

    def test_missing_S3_storage_credentials_2(self):
        expected_output = (
            """Invalid S3 config, missing access_key or secret_key in config file"""
        )
        if not sys.platform == "linux":
            pytest.skip(
                "for GHA, only works on linux (https://docs.github.com/en/actions/using-containerized-services/about-service-containers#about-service-containers)",
                allow_module_level=True,
            )
        with TmpMinIOStorage(minio_testcontainer) as tmp:
            storage_s3_json = benchmark_data_path / "storage_s3.json"
            auth_options = tmp.auth_options.copy()
            del auth_options["access_key"]
            print(auth_options)
            json.dump(auth_options, storage_s3_json.open("w"))

            os.environ["OB_STORAGE_S3_CONFIG"] = str(storage_s3_json)
            os.environ.pop("OB_STORAGE_S3_ACCESS_KEY", None)
            os.environ.pop("OB_STORAGE_S3_SECRET_KEY", None)

            _ = MinIOStorage(auth_options=tmp.auth_options, benchmark="benchmark004")

            subprocess.run(["ob", "run", "benchmark", "--benchmark", benchmark_path])
            with OmniCLISetup() as omni:
                result = omni.call(
                    [
                        "storage",
                        "create-version",
                        "--benchmark",
                        str(benchmark_path),
                    ]
                )
                assert result.exit_code == 1
                assert clean(result.output).startswith(clean(expected_output))


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
