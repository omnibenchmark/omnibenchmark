import sys

import pytest

import omni.cli.benchmark as ocb
from omni.io.MinIOStorage import MinIOStorage
from tests.io.MinIOStorage_setup import MinIOSetup, TmpMinIOStorage

# setup and start minio container
minio_testcontainer = MinIOSetup()


@pytest.mark.skipif(
    sys.platform != "linux", reason="for GHA, skip tests on non Linux platforms"
)
def test_list_benchmark(capfd):
    with TmpMinIOStorage(minio_testcontainer) as tmp:
        _ = MinIOStorage(auth_options=tmp.auth_options, benchmark="bm")
        ocb.list_benchmarks(tmp.auth_options_readonly["endpoint"])
        out, err = capfd.readouterr()
        assert len(out) > 0
