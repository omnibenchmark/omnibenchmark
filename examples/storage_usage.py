import json

from omnibenchmark.remote.RemoteStorage import StorageOptions
from omnibenchmark.remote.S3Storage import S3CompatibleStorage

with open("<CONFIG>.json", "r") as file:
    auth_options = json.load(file)

########################################################
### Example step-by-step usage of the S3CompatibleStorage class
########################################################

storage_options = StorageOptions(out_dir="out")

# Setup object — connects to an existing benchmark bucket or creates a new one
ms = S3CompatibleStorage(
    auth_options, benchmark="obob", storage_options=storage_options
)

# Upload objects to the benchmark bucket (boto3 keyword-argument style)
ms.client.put_object(Bucket=ms.benchmark, Key="out/f1.txt", Body=b"f1")
ms.client.put_object(Bucket=ms.benchmark, Key="out/f2.txt", Body=b"f2")

# Refresh the list of available benchmark versions
ms._get_versions()

# Set the version to create (None → auto-increment from latest)
ms.set_version()
print(f"Next version: {ms.version}")

# Or set a specific version manually
ms.set_version("0.1")
print(f"Version: {ms.version}")

# Create version 0.1 (tags current objects and writes the manifest)
ms.create_new_version()
print(f"Available versions: {ms.versions}")

# Update an object and create the next version
ms.client.put_object(Bucket=ms.benchmark, Key="out/f1.txt", Body=b"f1v2")
ms.set_version()  # auto-increments to 0.2
ms.create_new_version()

# Load and inspect objects for version 0.2
ms.load_objects()
print(ms.files)

# Load and inspect objects for version 0.1
ms.set_version("0.1")
ms.load_objects()
print(ms.files)

# Download an object to a local file
object_name = list(ms.files.keys())[0]
ms.download_object(object_name, "tmp.txt")
