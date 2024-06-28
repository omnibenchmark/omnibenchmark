import json

from omni.io.MinIOStorage import MinIOStorage

with open("<CONFIG>.json", "r") as file:
    auth_options = json.load(file)

########################################################
### Example step-by-step usage of the MinIOStorage class
########################################################

# setup object, existing benchmark or new benchmark
ms = MinIOStorage(auth_options, benchmark="tb")

# manually update overview bucket, usually not required
ms._update_overview()

# util: list containers
ms._get_containers()
# util: parse and list benchmarks
ms._get_benchmarks()
# util: get available benchmark versions
ms._get_versions()

# set version (latest)
ms.set_current_version()
print(f"{ms.major_version}.{ms.minor_version}")
# set version (manual)
ms.set_current_version(0, 1)
print(f"{ms.major_version}.{ms.minor_version}")

# parse new minor version (increment)
ms._parse_new_version(None, True)
# parse new major version (increment)
ms._parse_new_version(True, None)
# parse new version (manual)
ms._parse_new_version(0, 2)

# create new version
ms.set_new_version(0, 2)
# newly created version
print(f"{ms.major_version_new}.{ms.minor_version_new}")

# create new container for new version
ms._create_new_version()
# list versions of benchmark
ms._get_versions()

# list objects of version
ms._get_objects()
ms.files

# mark objects to copy to new version
ms.find_objects_to_copy()
ms.files

# copy flagged objects to new version
ms.copy_objects()
ms.files

########################################################
### Example one-step usage of the SwiftStorage class
########################################################

# setup object, existing benchmark or new benchmark
ms = MinIOStorage(auth_options, benchmark="tb")

ms.create_new_version()

########################################################
### Example public link download of object
########################################################

# auth only need public URL (key 'preauthurl')
auth_options_public = auth_options.copy()
del auth_options_public["access_key"]
del auth_options_public["secret_key"]
# retrieves versions from container tb2.overview which contains a list of versions
ms = MinIOStorage(auth_options_public, benchmark="tb")
ms.versions

# set version (latest)
ms.set_current_version()
ms.major_version
ms.minor_version

# list objects of version
ms._get_objects()
ms.files

# get test container content
ms.major_version = "test"
ms._get_objects()
ms.files

# get file content
names = list(ms.files.keys())
import requests

response = requests.get(
    f"https://{auth_options_public['endpoint']}/{ms.benchmark}.{ms.major_version}.{ms.minor_version}/{names[0]}"
)
response.text
