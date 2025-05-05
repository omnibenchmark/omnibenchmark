from itertools import groupby
from typing import Dict

import minio

from omnibenchmark.io.exception import MinIOStorageBucketManipulationException
from omnibenchmark.io.RemoteStorage import is_valid_version


def get_s3_object_versions_and_tags(
    client: minio.Minio, benchmark: str, readonly: bool = False
) -> Dict:
    """
    Retrieve the metadata of all objects in a S3 Bucket.

    Args:
    - client: A MinIO client object.
    - benchmark: The name of the S3 bucket.
    - readonly: A boolean indicating whether the client is read-only. Object tags can only be retrieved if the client is not read-only.

    Returns:
    - A dictionary of objects and associated metadata.
    """
    if not client.bucket_exists(benchmark):
        raise MinIOStorageBucketManipulationException(
            f"Benchmark {benchmark} does not exist."
        )

    object_ls = list(
        client.list_objects(
            benchmark, include_version=True, include_user_meta=True, recursive=True
        )
    )

    # Group listed objects by their name
    grouped_objects = groupby(object_ls, key=lambda o: o.object_name)

    di = {}

    # Iterate over object groups
    for object_name, objects in grouped_objects:
        objects = list(objects)
        di[object_name] = {}

        # Extract attributes for the objects group
        version_ids = [o.version_id for o in objects]
        sizes = [o.size for o in objects]
        last_modified = [o.last_modified for o in objects]
        is_delete_markers = [o.is_delete_marker for o in objects]
        etags = [o.etag for o in objects]

        # Iterate to all versions and create a dictionary entry
        for i, v in enumerate(version_ids):
            tags_filt = {}
            if not readonly and not is_delete_markers[i]:
                tags = client.get_object_tags(benchmark, object_name, version_id=v)
                if tags is not None:
                    tags_filt = {k: w for k, w in tags.items() if is_valid_version(k)}

            di[object_name][v] = {
                "tags": tags_filt,
                "size": sizes[i],
                "last_modified": last_modified[i],
                "is_delete_marker": is_delete_markers[i],
                "etag": etags[i],
            }
    return di
