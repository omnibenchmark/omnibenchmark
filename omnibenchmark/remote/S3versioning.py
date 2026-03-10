"""S3 object versioning and tagging helpers (boto3-based)."""

from __future__ import annotations

from typing import Dict

try:
    from botocore.exceptions import ClientError
except ImportError:
    ClientError = Exception  # type: ignore[assignment,misc]

from omnibenchmark.remote.exception import MinIOStorageBucketManipulationException
from omnibenchmark.remote.RemoteStorage import is_valid_version


def get_s3_object_versions_and_tags(
    client, benchmark: str, readonly: bool = False
) -> Dict:
    """
    Retrieve the metadata of all objects in an S3 bucket.

    Args:
        client: A boto3 S3 client.
        benchmark: The name of the S3 bucket.
        readonly: When True, skip fetching object tags (reduces API calls).

    Returns:
        A nested dict keyed by object name → version_id → metadata dict with
        keys: tags, size, last_modified, is_delete_marker, etag.
    """
    try:
        client.head_bucket(Bucket=benchmark)
    except ClientError:
        raise MinIOStorageBucketManipulationException(
            f"Benchmark {benchmark} does not exist."
        )

    paginator = client.get_paginator("list_object_versions")
    di: Dict = {}

    try:
        for page in paginator.paginate(Bucket=benchmark):
            for version in page.get("Versions", []):
                key = version["Key"]
                vid = version["VersionId"]
                if key not in di:
                    di[key] = {}

                tags_filt: dict = {}
                if not readonly:
                    try:
                        tag_response = client.get_object_tagging(
                            Bucket=benchmark, Key=key, VersionId=vid
                        )
                        tags_filt = {
                            t["Key"]: t["Value"]
                            for t in tag_response.get("TagSet", [])
                            if is_valid_version(t["Key"])
                        }
                    except ClientError:
                        # Tag retrieval failures are non-fatal; keep tags empty and continue.
                        tags_filt = {}

                di[key][vid] = {
                    "tags": tags_filt,
                    "size": version["Size"],
                    "last_modified": version["LastModified"],
                    "is_delete_marker": False,
                    "etag": version["ETag"].strip('"'),
                }

            for dm in page.get("DeleteMarkers", []):
                key = dm["Key"]
                vid = dm["VersionId"]
                if key not in di:
                    di[key] = {}
                di[key][vid] = {
                    "tags": {},
                    "size": 0,
                    "last_modified": dm["LastModified"],
                    "is_delete_marker": True,
                    "etag": "",
                }

    except ClientError as e:
        code = getattr(e, "response", {}).get("Error", {}).get("Code", "")
        if code in ("NoSuchBucket", "404"):
            raise MinIOStorageBucketManipulationException(
                f"Benchmark {benchmark} does not exist."
            )
        raise

    return di
