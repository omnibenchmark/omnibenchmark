import pytest
import requests
import os
import omni.io
import omni.io.utils as oiu


def test_get_storage():
    ss = oiu.get_storage(
        storage_type="s3",
        auth_options={
            "endpoint": "https://omnibenchmark.mls.uzh.ch:9000",
            "secure": False,
        },
        benchmark="bm",
    )
    assert isinstance(ss, omni.io.S3Storage.S3Storage)

    with pytest.raises(requests.exceptions.HTTPError):
        ss = oiu.get_storage(
            storage_type="s3",
            auth_options={
                "endpoint": "https://omnibenchmark.mls.uzh.ch:9000",
                "secure": False,
            },
            benchmark="not_existing_benchmark",
        )

    ss = oiu.get_storage(
        storage_type="minio",
        auth_options={
            "endpoint": "https://omnibenchmark.mls.uzh.ch:9000",
            "secure": False,
        },
        benchmark="bm",
    )
    assert isinstance(ss, omni.io.MinIOStorage.MinIOStorage)

    with pytest.raises(requests.exceptions.HTTPError):
        ss = oiu.get_storage(
            storage_type="minio",
            auth_options={
                "endpoint": "https://omnibenchmark.mls.uzh.ch:9000",
                "secure": False,
            },
            benchmark="not_existing_benchmark",
        )


def cleanup_md5():
    if os.path.exists("tests/io/md5sum_example.txt"):
        os.remove("tests/io/md5sum_example.txt")


@pytest.fixture(scope="session", autouse=True)
def cleanup(request):
    """Remove the file created during the test."""
    request.addfinalizer(cleanup_md5)


def test_md5():
    with open("tests/io/md5sum_example.txt", "w") as f:
        f.write("asdfasdf")
    assert oiu.md5("tests/io/md5sum_example.txt") == "6a204bd89f3c8348afd5c77c717a097a"

    with pytest.raises(FileNotFoundError):
        oiu.md5("not_existing_file.txt")


def test_sizeof_fmt():
    assert oiu.sizeof_fmt(0) == "    0B"
    assert oiu.sizeof_fmt(1) == "    1B"
    assert oiu.sizeof_fmt(1023) == " 1023B"
    assert oiu.sizeof_fmt(1024) == "1.0KiB"
    assert oiu.sizeof_fmt(1024**2) == "1.0MiB"
    assert oiu.sizeof_fmt(1024**3) == "1.0GiB"
    assert oiu.sizeof_fmt(1024**4) == "1.0TiB"
    assert oiu.sizeof_fmt(1024**5) == "1.0PiB"
    assert oiu.sizeof_fmt(1024**6) == "1.0EiB"
    assert oiu.sizeof_fmt(1024**7) == "1.0ZiB"
    assert oiu.sizeof_fmt(1024**8) == "1.0YiB"
