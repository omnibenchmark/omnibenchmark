import minio
import pytest

import omnibenchmark.io.utils as oiu
from omnibenchmark.io.RemoteStorage import RemoteStorage

from ..fixtures import minio_storage, _minio_container  # noqa: F401


# fmt: off
def test_get_storage_raise_exception_when_passed_invalid_benchmark_path(minio_storage):  # noqa: F811
# fmt: on
    # happy path
    storage = oiu.get_storage(
        storage_type="minio",
        auth_options=minio_storage.auth_options,
        benchmark="test",
    )
    assert isinstance(storage, RemoteStorage)

    # sad path
    with pytest.raises(minio.error.S3Error):
        oiu.get_storage(
            storage_type="minio",
            auth_options=minio_storage.auth_options,
            benchmark="not_existing_benchmark",
        )


# TODO: we should use humanize for this, it's already in the deps
@pytest.mark.short
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
