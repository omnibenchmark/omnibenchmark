import datetime
import io
import sys
from packaging.version import Version
from pathlib import Path

import pytest
import yaml

from omni_schema.datamodel.omni_schema import SoftwareBackendEnum
from omnibenchmark.benchmark import Benchmark
from omnibenchmark.io.exception import RemoteStorageInvalidInputException
from omnibenchmark.io.MinIOStorage import MinIOStorage
from tests.io.MinIOStorage_setup import MinIOSetup, TmpMinIOStorage

if not sys.platform == "linux":
    pytest.skip(
        "for GHA, only works on linux (https://docs.github.com/en/actions/using-containerized-services/about-service-containers#about-service-containers)",
        allow_module_level=True,
    )

# setup and start minio container
# TODO: no please, use session fixture
minio_testcontainer = MinIOSetup()  # sys.platform == "linux"


benchmark_data = Path("..") / "data"
benchmark_data_path = Path(__file__).parent / benchmark_data


class TestMinIOStorage:
    def test_init_fail(self):
        with pytest.raises(AssertionError):
            MinIOStorage(auth_options={}, benchmark="test")

    def test_init_success(self):
        with TmpMinIOStorage(minio_testcontainer) as tmp:
            _ = MinIOStorage(auth_options=tmp.auth_options, benchmark=tmp.bucket_name)
            _ = tmp.minio.get_client()
            ss = MinIOStorage(auth_options=tmp.auth_options, benchmark=tmp.bucket_name)
            assert ss.benchmark == tmp.bucket_name

    def test_init_public_success(self):
        with TmpMinIOStorage(minio_testcontainer) as tmp:
            _ = MinIOStorage(auth_options=tmp.auth_options, benchmark=tmp.bucket_name)
            ss = MinIOStorage(
                auth_options=tmp.auth_options_readonly, benchmark=tmp.bucket_name
            )
            assert ss.benchmark == tmp.bucket_name

    def test__test_connect_success_with_valid_endpoint(self):
        with TmpMinIOStorage(minio_testcontainer) as tmp:
            ss = MinIOStorage(auth_options=tmp.auth_options, benchmark=tmp.bucket_name)
            ss._test_connect()

    def test__create_benchmark_success_create_benchmark(self):
        with TmpMinIOStorage(minio_testcontainer) as tmp:
            ss = MinIOStorage(auth_options=tmp.auth_options, benchmark=tmp.bucket_name)
            ss._create_benchmark(f"{tmp.bucket_name}2")

            assert tmp.minio.get_client().bucket_exists(f"{tmp.bucket_name}")
            assert tmp.minio.get_client().bucket_exists(f"{tmp.bucket_name}2")

    def test__get_versions_success_get_version(self):
        with TmpMinIOStorage(minio_testcontainer) as tmp:
            ss = MinIOStorage(auth_options=tmp.auth_options, benchmark=tmp.bucket_name)
            ss.set_version("0.1")
            ss.create_new_version()
            ss._get_versions()
            ss = MinIOStorage(auth_options=tmp.auth_options, benchmark=tmp.bucket_name)
            assert ss.versions == [Version("0.1")]

    def test__get_versions_public_success_get_version(self):
        with TmpMinIOStorage(minio_testcontainer) as tmp:
            ss = MinIOStorage(auth_options=tmp.auth_options, benchmark=tmp.bucket_name)
            ss.set_version("0.1")
            ss.create_new_version()
            ss = MinIOStorage(
                auth_options=tmp.auth_options_readonly, benchmark=tmp.bucket_name
            )
            ss._get_versions()
            assert ss.versions == [Version("0.1")]

    def test__create_new_version_fail_if_no_version_provided(self):
        with TmpMinIOStorage(minio_testcontainer) as tmp:
            ss = MinIOStorage(auth_options=tmp.auth_options, benchmark=tmp.bucket_name)
        with pytest.raises(RemoteStorageInvalidInputException):
            ss.create_new_version()

    def test__create_new_version_success_if_version_provided(self):
        with TmpMinIOStorage(minio_testcontainer) as tmp:
            ss = MinIOStorage(auth_options=tmp.auth_options, benchmark=tmp.bucket_name)
            ss.set_version("0.1")
            ss.create_new_version()
            ss.set_version()
            ss.create_new_version()
            assert Version("0.2") in ss.versions

    def test__create_new_version_success_if_version_provided2(self):
        with TmpMinIOStorage(minio_testcontainer) as tmp:
            MinIOStorage(auth_options=tmp.auth_options, benchmark=tmp.bucket_name)

    def test__get_objects(self):
        with TmpMinIOStorage(minio_testcontainer) as tmp:
            ss = MinIOStorage(auth_options=tmp.auth_options, benchmark=tmp.bucket_name)
            ss.client.put_object(ss.benchmark, "out/file1.txt", io.BytesIO(b""), 0)
            ss.client.put_object(ss.benchmark, "out/file2.txt", io.BytesIO(b""), 0)
            ss.set_version()
            ss._get_objects()
            assert all(
                [f in ss.files.keys() for f in ["out/file1.txt", "out/file2.txt"]]
            )
            assert ss.files["out/file1.txt"].keys() == {
                "version_id",
                "etag",
                "last_modified",
                "size",
            }
            assert ss.files["out/file2.txt"].keys() == {
                "version_id",
                "etag",
                "last_modified",
                "size",
            }
            assert ss.files["out/file1.txt"]["size"] == 0
            assert ss.files["out/file2.txt"]["size"] == 0
            assert (
                ss.files["out/file1.txt"]["etag"] == "d41d8cd98f00b204e9800998ecf8427e"
            )
            assert (
                ss.files["out/file2.txt"]["etag"] == "d41d8cd98f00b204e9800998ecf8427e"
            )
            assert isinstance(
                ss.files["out/file1.txt"]["last_modified"], datetime.datetime
            )
            assert isinstance(
                ss.files["out/file2.txt"]["last_modified"], datetime.datetime
            )

    def test__get_objects_public(self):
        with TmpMinIOStorage(minio_testcontainer) as tmp:
            ss = MinIOStorage(auth_options=tmp.auth_options, benchmark=tmp.bucket_name)
            ss.client.put_object(ss.benchmark, "out/file1.txt", io.BytesIO(b""), 0)
            ss.client.put_object(ss.benchmark, "out/file2.txt", io.BytesIO(b""), 0)
            ss = MinIOStorage(
                auth_options=tmp.auth_options_readonly, benchmark=tmp.bucket_name
            )
            ss.set_version()
            ss._get_objects()
            assert all(
                [f in ss.files.keys() for f in ["out/file1.txt", "out/file2.txt"]]
            )
            assert ss.files["out/file1.txt"].keys() == {
                "version_id",
                "etag",
                "last_modified",
                "size",
            }
            assert ss.files["out/file2.txt"].keys() == {
                "version_id",
                "etag",
                "last_modified",
                "size",
            }
            assert ss.files["out/file1.txt"]["size"] == 0
            assert ss.files["out/file2.txt"]["size"] == 0
            assert (
                ss.files["out/file1.txt"]["etag"] == "d41d8cd98f00b204e9800998ecf8427e"
            )
            assert (
                ss.files["out/file2.txt"]["etag"] == "d41d8cd98f00b204e9800998ecf8427e"
            )
            assert isinstance(
                ss.files["out/file1.txt"]["last_modified"], datetime.datetime
            )
            assert isinstance(
                ss.files["out/file2.txt"]["last_modified"], datetime.datetime
            )

    def test_create_new_version(self):
        with TmpMinIOStorage(minio_testcontainer) as tmp:
            ss = MinIOStorage(auth_options=tmp.auth_options, benchmark=tmp.bucket_name)
            _ = ss.client.put_object(ss.benchmark, "out/file1.txt", io.BytesIO(b""), 0)
            _ = ss.client.put_object(ss.benchmark, "out/file2.txt", io.BytesIO(b""), 0)
            ss.set_version("0.2")
            ss.create_new_version()
            ss._get_objects()
            assert all(
                [f in ss.files.keys() for f in ["out/file1.txt", "out/file2.txt"]]
            )
            assert ss.files["out/file1.txt"].keys() == {
                "version_id",
                "etag",
                "last_modified",
                "size",
            }
            assert ss.files["out/file2.txt"].keys() == {
                "version_id",
                "etag",
                "last_modified",
                "size",
            }

    def test_filter_with_benchmark(self):
        with TmpMinIOStorage(minio_testcontainer) as tmp:
            ss = MinIOStorage(auth_options=tmp.auth_options, benchmark=tmp.bucket_name)
            _ = ss.client.put_object(ss.benchmark, "out/file1.txt", io.BytesIO(b""), 0)
            _ = ss.client.put_object(ss.benchmark, "out/file2.txt", io.BytesIO(b""), 0)
            ss.set_version("0.2")
            ss.create_new_version()
            ss._get_objects()
            ss.files
            benchmark = str(benchmark_data_path / "mock_benchmark.yaml")
            print(benchmark)
            with open(benchmark, "r") as fh:
                yaml.safe_load(fh)
                benchmark = Benchmark(Path(benchmark))
            ss.set_version("0.3")
            ss.create_new_version(benchmark)
            ss._get_objects()
            assert all(
                [f not in ss.files.keys() for f in ["out/file1.txt", "out/file2.txt"]]
            )

    def test_store_software_and_config_with_benchmark(self):
        with TmpMinIOStorage(minio_testcontainer) as tmp:
            ss = MinIOStorage(auth_options=tmp.auth_options, benchmark=tmp.bucket_name)
            benchmark = str(benchmark_data_path / "Clustering.yaml")
            print(benchmark)
            with open(benchmark, "r") as fh:
                yaml.safe_load(fh)
                benchmark = Benchmark(Path(benchmark))

            benchmark.converter.model.software_backend = SoftwareBackendEnum("conda")
            ss.set_version("0.3")
            ss.create_new_version(benchmark)
            ss._get_objects()
            ss.files.keys()
            assert all(
                [
                    f in ss.files.keys()
                    for f in [
                        "software/R_4.4.1_Clustering.yaml",
                        "software/Python_3.12.6_Clustering.yaml",
                        "config/benchmark.yaml",
                    ]
                ]
            )


def cleanup_buckets_on_exit():
    """Cleanup a testing directory once we are finished."""
    TmpMinIOStorage(minio_testcontainer).cleanup_buckets()


@pytest.fixture(scope="session", autouse=True)
def cleanup(request):
    """Cleanup a testing directory once we are finished."""
    request.addfinalizer(cleanup_buckets_on_exit)
