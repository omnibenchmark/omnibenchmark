import omni.io.files as oif
import json
import requests
import os
from pathlib import Path
import pytest

if Path("tests/.minio_test_config.json").exists():
    with open("tests/.minio_test_config.json", "r") as file:
        auth_options = json.load(file)
elif "S3_ENDPOINT_URL" in os.environ:
    auth_options = {}
    auth_options["endpoint"] = os.environ.get("S3_ENDPOINT_URL")
    auth_options["secure"] = False
else:
    raise ValueError("No S3 credentials found")

endpoint = (
    f"http://{auth_options['endpoint'].replace('http://', '').replace('https://', '')}"
)


def test_get_benchmarks_public():
    benchmarks = oif.get_benchmarks_public(endpoint)
    assert isinstance(benchmarks, list)

    with pytest.raises(requests.exceptions.MissingSchema):
        benchmarks = oif.get_benchmarks_public("non_existing_endpoint")

    with pytest.raises(requests.exceptions.ConnectionError):
        benchmarks = oif.get_benchmarks_public("http://non_existing_endpoint")


def test_get_benchmark_versions_public():
    versions = oif.get_benchmark_versions_public("bm", endpoint)

    assert isinstance(versions, list)

    with pytest.raises(requests.exceptions.MissingSchema):
        versions = oif.get_benchmark_versions_public("bm", "non_existing_endpoint")

    with pytest.raises(requests.exceptions.ConnectionError):
        versions = oif.get_benchmark_versions_public(
            "bm", "http://non_existing_endpoint"
        )

    with pytest.raises(requests.exceptions.HTTPError):
        versions = oif.get_benchmark_versions_public("non_existing_benchmark", endpoint)
