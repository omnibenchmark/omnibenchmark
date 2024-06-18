import sys

import pytest
import requests

import omni.io.files as oif
from omni.io.MinIOStorage import MinIOStorage
from tests.io.MinIOStorage_setup import MinIOSetup, TmpMinIOStorage

# setup and start minio container
minio_testcontainer = MinIOSetup()


@pytest.mark.skipif(
    sys.platform != "linux", reason="for GHA, skip tests on non Linux platforms"
)
def test_get_benchmarks_public():
    with TmpMinIOStorage(minio_testcontainer) as tmp:
        _ = MinIOStorage(auth_options=tmp.auth_options, benchmark="bm")
        benchmarks = oif.get_benchmarks_public(tmp.auth_options_readonly["endpoint"])
        assert isinstance(benchmarks, list)

        with pytest.raises(requests.exceptions.MissingSchema):
            benchmarks = oif.get_benchmarks_public("non_existing_endpoint")

        with pytest.raises(requests.exceptions.ConnectionError):
            benchmarks = oif.get_benchmarks_public("http://non_existing_endpoint")


@pytest.mark.skipif(
    sys.platform != "linux", reason="for GHA, skip tests on non Linux platforms"
)
def test_get_benchmark_versions_public():
    with TmpMinIOStorage(minio_testcontainer) as tmp:
        _ = MinIOStorage(auth_options=tmp.auth_options, benchmark="bm")
        versions = oif.get_benchmark_versions_public(
            "bm", tmp.auth_options_readonly["endpoint"]
        )

        assert isinstance(versions, list)

        with pytest.raises(requests.exceptions.MissingSchema):
            versions = oif.get_benchmark_versions_public("bm", "non_existing_endpoint")

        with pytest.raises(requests.exceptions.ConnectionError):
            versions = oif.get_benchmark_versions_public(
                "bm", "http://non_existing_endpoint"
            )

        with pytest.raises(requests.exceptions.HTTPError):
            versions = oif.get_benchmark_versions_public(
                "non_existing_benchmark", tmp.auth_options_readonly["endpoint"]
            )
