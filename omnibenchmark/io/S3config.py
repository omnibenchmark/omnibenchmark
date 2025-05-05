import json
import os
import sys


from omnibenchmark.cli.utils.logging import logger
from omnibenchmark.io.RemoteStorage import DEFAULT_STORAGE_OPTIONS


class S3_DEFAULT_STORAGE_OPTIONS(DEFAULT_STORAGE_OPTIONS):
    pass


def benchmarker_access_token_policy(benchmark):
    """S3 policy for access token for specific benchmark, allows archiving (Governance Retention)"""
    return {
        "Version": "2012-10-17",
        "Statement": [
            {
                "Effect": "Allow",
                "Action": ["s3:*"],
                "Resource": [
                    f"arn:aws:s3:::{benchmark}/*",
                    f"arn:aws:s3:::{benchmark}",
                ],
            },
            {
                "Effect": "Deny",
                "Action": [
                    "s3:BypassGovernanceRetention",
                    "s3:DeleteObjectTagging",
                    "s3:DeleteObjectVersion",
                    "s3:DeleteObjectVersionTagging",
                    "s3:DeleteBucket",
                    "s3:DeleteBucketPolicy",
                ],
                "Resource": [
                    f"arn:aws:s3:::{benchmark}/*",
                    f"arn:aws:s3:::{benchmark}",
                ],
            },
        ],
    }


def bucket_readonly_policy(bucket_name):
    return {
        "Version": "2012-10-17",
        "Statement": [
            {
                "Effect": "Allow",
                "Principal": {"AWS": "*"},
                "Action": ["s3:GetBucketLocation", "s3:ListBucket"],
                "Resource": f"arn:aws:s3:::{bucket_name}",
            },
            {
                "Effect": "Allow",
                "Principal": {"AWS": "*"},
                "Action": "s3:GetObject",
                "Resource": f"arn:aws:s3:::{bucket_name}/*",
            },
        ],
    }


def S3_access_config_from_env() -> dict:
    """Get S3 access config from environment variables or file"""
    if (
        "OB_STORAGE_S3_ACCESS_KEY" in os.environ
        and "OB_STORAGE_S3_SECRET_KEY" in os.environ
    ):
        return {
            "access_key": os.environ["OB_STORAGE_S3_ACCESS_KEY"],
            "secret_key": os.environ["OB_STORAGE_S3_SECRET_KEY"],
        }
    elif "OB_STORAGE_S3_CONFIG" in os.environ:
        with open(os.environ["OB_STORAGE_S3_CONFIG"], "r") as file:
            auth_options = json.load(file)
        if "access_key" in auth_options and "secret_key" in auth_options:
            return auth_options
        else:
            logger.error(
                f"Invalid S3 config, missing access_key or secret_key in config file ({os.environ['OB_STORAGE_S3_CONFIG']})",
            )
            sys.exit(1)
    else:
        logger.error(
            "Invalid S3 config. Missing access_key and secret_key in environment variables (OB_STORAGE_S3_ACCESS_KEY, OB_STORAGE_S3_SECRET_KEY) or OB_STORAGE_S3_CONFIG",
        )
        sys.exit(1)


if __name__ == "__main__":
    print(json.dumps(benchmarker_access_token_policy("obob"), indent=2))
    print(json.dumps(bucket_readonly_policy("omnibenchmarktestbucket"), indent=2))
