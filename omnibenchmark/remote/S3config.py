import json
import os
import sys
from pathlib import Path

from dotenv import load_dotenv

from omnibenchmark.cli.utils.logging import logger


# Load .env file from the current working directory if it exists
# This allows users to store S3 credentials in a .env file
dotenv_path = Path.cwd() / ".env"
if dotenv_path.exists():
    load_dotenv(dotenv_path=dotenv_path, override=False)
else:
    # Try loading from the project root (in case we're running from a subdirectory)
    # Walk up the directory tree to find .env
    current_path = Path.cwd()
    for parent in [current_path] + list(current_path.parents):
        dotenv_path = parent / ".env"
        if dotenv_path.exists():
            load_dotenv(dotenv_path=dotenv_path, override=False)
            break


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


def S3_access_config_from_env(required: bool = True) -> dict:
    """Get S3 access config from environment variables or file

    Args:
        required: If True, raise error when credentials are missing.
                 If False, return empty dict when credentials are missing.

    Returns:
        dict: Dictionary with access_key and secret_key, or empty dict if not required
    """
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
            if required:
                import click

                logger.error(
                    click.style("[ERROR]", fg="red", bold=True)
                    + f" Missing access_key or secret_key in config file: {os.environ['OB_STORAGE_S3_CONFIG']}"
                )
                sys.exit(1)
            return {}
    else:
        if required:
            import click

            logger.error(
                click.style("[ERROR]", fg="red", bold=True)
                + " Missing S3 credentials. Set OB_STORAGE_S3_ACCESS_KEY and OB_STORAGE_S3_SECRET_KEY, or OB_STORAGE_S3_CONFIG"
            )
            sys.exit(1)
        return {}


if __name__ == "__main__":
    print(json.dumps(benchmarker_access_token_policy("obob"), indent=2))
    print(json.dumps(bucket_readonly_policy("omnibenchmarktestbucket"), indent=2))
