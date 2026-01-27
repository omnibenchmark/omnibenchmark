import datetime
import io
import tempfile
import shutil
from packaging.version import Version
from pathlib import Path

import pytest
from git import Repo

from omnibenchmark.model import SoftwareBackendEnum
from omnibenchmark.benchmark import BenchmarkExecution
from omnibenchmark.remote.RemoteStorage import StorageOptions
from omnibenchmark.remote.exception import RemoteStorageInvalidInputException
from omnibenchmark.remote.MinIOStorage import MinIOStorage

from ..fixtures import minio_storage, _minio_container  # noqa: F401


def get_benchmark_data_path() -> Path:
    return Path(__file__).resolve().parent.parent / "data"


class TestMinIOStorage:
    def test_init_fail(self):
        with pytest.raises(AssertionError):
            MinIOStorage(
                auth_options={},
                benchmark="test",
                storage_options=StorageOptions(out_dir="out"),
            )

    def test_init_success(self, minio_storage):  # noqa: F811
        client = minio_storage.get_storage_client()
        assert client.benchmark == minio_storage.bucket_name

    # fmt: off
    def test__test_connect_success_with_valid_endpoint(self, minio_storage):  # noqa: F811
    # fmt: on
        client = minio_storage.get_storage_client()
        client._test_connect()

    # fmt: off
    def test__create_benchmark_success_create_benchmark(self, minio_storage):  # noqa: F811
    # fmt: on
        client = minio_storage.get_storage_client()
        # XXX should not test private methods
        client._create_benchmark(f"{minio_storage.bucket_name}2")

        assert client.client.bucket_exists(f"{minio_storage.bucket_name}")
        assert client.client.bucket_exists(f"{minio_storage.bucket_name}2")

    # fmt: off
    def test__get_versions_success_get_version(self, minio_storage):  # noqa: F811
    # fmt: on
        client = minio_storage.get_storage_client()
        client.set_version("0.1")
        client.create_new_version()
        client2 = minio_storage.get_storage_client()
        assert client2.versions == [Version("0.1")]

    # fmt: off
    def test__get_versions_public_success_get_version(self, minio_storage):  # noqa: F811
    # fmt: on
        client = minio_storage.get_storage_client()
        client.set_version("0.1")
        client.create_new_version()

        client2 = minio_storage.get_storage_client()
        assert client2.versions == [Version("0.1")]

    # fmt: off
    def test__create_new_version_fail_if_no_version_provided(self, minio_storage):  # noqa: F811
    # fmt: on
        client = minio_storage.get_storage_client()
        with pytest.raises(RemoteStorageInvalidInputException):
            client.create_new_version()

    # fmt: off
    def test__create_new_version_success_if_version_provided(self, minio_storage):  # noqa: F811
    # fmt: on
        client = minio_storage.get_storage_client()
        client.set_version("0.1")
        client.create_new_version()
        client.set_version()
        client.create_new_version()
        assert Version("0.1") in client.versions
        assert Version("0.2") in client.versions

    # fmt: off
    def test__get_objects(self, minio_storage):  # noqa: F811
    # fmt: on
        client = minio_storage.get_storage_client()
        client.client.put_object(client.benchmark, "out/file1.txt", io.BytesIO(b""), 0)
        client.client.put_object(client.benchmark, "out/file2.txt", io.BytesIO(b""), 0)
        client.set_version()
        client._get_objects()
        assert all(
            [f in client.files.keys() for f in ["out/file1.txt", "out/file2.txt"]]
        )

        assert client.files["out/file1.txt"].keys() == {
            "version_id",
            "etag",
            "last_modified",
            "size",
        }
        assert client.files["out/file2.txt"].keys() == {
            "version_id",
            "etag",
            "last_modified",
            "size",
        }
        assert client.files["out/file1.txt"]["size"] == 0
        assert client.files["out/file2.txt"]["size"] == 0
        assert (
            client.files["out/file1.txt"]["etag"] == "d41d8cd98f00b204e9800998ecf8427e"
        )
        assert (
            client.files["out/file2.txt"]["etag"] == "d41d8cd98f00b204e9800998ecf8427e"
        )
        assert isinstance(
            client.files["out/file1.txt"]["last_modified"], datetime.datetime
        )
        assert isinstance(
            client.files["out/file2.txt"]["last_modified"], datetime.datetime
        )

    # fmt: off
    def test__get_objects_public(self, minio_storage):  # noqa: F811
    # fmt: on
        client = minio_storage.get_storage_client()
        client.client.put_object(client.benchmark, "out/file1.txt", io.BytesIO(b""), 0)
        client.client.put_object(client.benchmark, "out/file2.txt", io.BytesIO(b""), 0)
        client.set_version()
        client._get_objects()
        assert all(
            [f in client.files.keys() for f in ["out/file1.txt", "out/file2.txt"]]
        )
        assert client.files["out/file1.txt"].keys() == {
            "version_id",
            "etag",
            "last_modified",
            "size",
        }
        assert client.files["out/file2.txt"].keys() == {
            "version_id",
            "etag",
            "last_modified",
            "size",
        }
        assert client.files["out/file1.txt"]["size"] == 0
        assert client.files["out/file2.txt"]["size"] == 0
        assert (
            client.files["out/file1.txt"]["etag"] == "d41d8cd98f00b204e9800998ecf8427e"
        )
        assert (
            client.files["out/file2.txt"]["etag"] == "d41d8cd98f00b204e9800998ecf8427e"
        )
        assert isinstance(
            client.files["out/file1.txt"]["last_modified"], datetime.datetime
        )
        assert isinstance(
            client.files["out/file2.txt"]["last_modified"], datetime.datetime
        )

    # fmt: off
    def test_create_new_version(self, minio_storage):  # noqa: F811
    # fmt: on
        client = minio_storage.get_storage_client()
        _ = client.client.put_object(
            client.benchmark, "out/file1.txt", io.BytesIO(b""), 0
        )
        _ = client.client.put_object(
            client.benchmark, "out/file2.txt", io.BytesIO(b""), 0
        )
        client.set_version("0.2")
        client.create_new_version()
        client._get_objects()
        assert all(
            [f in client.files.keys() for f in ["out/file1.txt", "out/file2.txt"]]
        )
        assert client.files["out/file1.txt"].keys() == {
            "version_id",
            "etag",
            "last_modified",
            "size",
        }
        assert client.files["out/file2.txt"].keys() == {
            "version_id",
            "etag",
            "last_modified",
            "size",
        }

    # fmt: off
    def test_filter_with_benchmark(self, minio_storage):  # noqa: F811
    # fmt: on
        client = minio_storage.get_storage_client()
        _ = client.client.put_object(
            client.benchmark, "out/file1.txt", io.BytesIO(b""), 0
        )
        _ = client.client.put_object(
            client.benchmark, "out/file2.txt", io.BytesIO(b""), 0
        )
        client.set_version("0.2")
        client.create_new_version()
        client._get_objects()
        client.files

        # Set up a temporary git repository for testing version persistence
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            # Initialize git repo
            repo = Repo.init(temp_path)

            # Copy benchmark file to git repo
            path = get_benchmark_data_path()
            original_benchmark_file = path / "mock_benchmark.yaml"
            benchmark_file = temp_path / "mock_benchmark.yaml"
            shutil.copy2(original_benchmark_file, benchmark_file)

            # Commit the initial benchmark file
            repo.index.add([str(benchmark_file.name)])
            repo.index.commit("Initial benchmark commit")

            # Create BenchmarkExecution with git-tracked file
            benchmark = BenchmarkExecution(benchmark_file)

            client.set_version("0.3")
            client.create_new_version(benchmark)
            client._get_objects()
            assert all(
                [f not in client.files.keys() for f in ["out/file1.txt", "out/file2.txt"]]
            )

    # fmt: off
    def test_store_software_and_config_with_benchmark(self, minio_storage):  # noqa: F811
    # fmt: on

        client = minio_storage.get_storage_client()

        # Set up a temporary git repository for testing version persistence
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            # Initialize git repo
            repo = Repo.init(temp_path)

            # Copy benchmark file to git repo
            path = get_benchmark_data_path()
            original_benchmark_file = path / "Clustering.yaml"
            benchmark_file = temp_path / "Clustering.yaml"
            shutil.copy2(original_benchmark_file, benchmark_file)

            # Copy environment files as well
            envs_dir = temp_path / "envs"
            envs_dir.mkdir()
            original_envs_dir = path / "envs"
            for env_file in original_envs_dir.glob("*.yaml"):
                shutil.copy2(env_file, envs_dir / env_file.name)

            # Commit the initial files
            repo.index.add([str(benchmark_file.name), "envs/*"])
            repo.index.commit("Initial benchmark commit")

            # Create BenchmarkExecution with git-tracked file
            benchmark = BenchmarkExecution(benchmark_file)

            benchmark.model.software_backend = SoftwareBackendEnum("conda")

            client.set_version("0.3")
            client.create_new_version(benchmark)
            client._get_objects()
            client.files.keys()
            assert all(
                [
                    f in client.files.keys()
                    for f in [
                        "software/R_4.4.1_Clustering.yaml",
                        "software/Python_3.12.6_Clustering.yaml",
                        "config/benchmark.yaml",
                    ]
                ]
            )
