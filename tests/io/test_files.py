import sys

from tests.io.MinIOStorage_setup import MinIOSetup

# setup and start minio container
minio_testcontainer = MinIOSetup(sys.platform == "linux")
