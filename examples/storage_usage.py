import io
import json

from omni.io.MinIOStorage import MinIOStorage

with open("<CONFIG>.json", "r") as file:
    auth_options = json.load(file)

########################################################
### Example step-by-step usage of the MinIOStorage class
########################################################

# setup object, existing benchmark or new benchmark
ms = MinIOStorage(auth_options, benchmark="obob")

# upload objects to benchmark
_ = ms.client.put_object(ms.benchmark, "out/f1.txt", io.BytesIO(b"f1"), 2)
_ = ms.client.put_object(ms.benchmark, "out/f2.txt", io.BytesIO(b"f2"), 2)

# util: get available benchmark versions
ms._get_versions()

# set version (latest)
ms.set_version()
# version to create
print(f"{ms.version}")
# set version (manual)
ms.set_version("0.1")
# version to create
print(f"{ms.version}")

# create version
ms.create_new_version()
# list versions of benchmark
ms.versions

# update object
_ = ms.client.put_object(ms.benchmark, "out/f1.txt", io.BytesIO(b"f1v2"), 4)
# prepare new version
ms.set_version()
# create version
ms.create_new_version()

# list objects of version 0.2
ms._get_objects()
ms.files

# list objects of version 0.1
ms.set_version("0.1")
ms._get_objects()
ms.files


# download object to file (tmp.txt)
object_name = list(ms.files.keys())[0]
ms.download_object(object_name, "tmp.txt")
