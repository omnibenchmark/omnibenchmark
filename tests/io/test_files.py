import sys

import pytest
import requests

import omni.io.files as oif
from omni.io.MinIOStorage import MinIOStorage
from tests.io.MinIOStorage_setup import MinIOSetup, TmpMinIOStorage

# setup and start minio container
minio_testcontainer = MinIOSetup(sys.platform == "linux")
