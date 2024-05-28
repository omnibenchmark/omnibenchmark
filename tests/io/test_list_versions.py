import omni.cli.benchmark as ocb
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


def test_list_benchmark(capfd):
    ocb.list_versions("bm", endpoint)
    out, err = capfd.readouterr()
    assert len(out) > 0

    with pytest.raises(requests.exceptions.HTTPError):
        ocb.list_versions("not_existing_benchmark", endpoint)
