from omni.io.RemoteStorage import DEFAULT_STORAGE_OPTIONS

S3_DEFAULT_STORAGE_OPTIONS = {**DEFAULT_STORAGE_OPTIONS}


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


if __name__ == "__main__":
    import json

    print(json.dumps(benchmarker_access_token_policy("obob"), indent=2))
    print(json.dumps(bucket_readonly_policy("omnibenchmarktestbucket"), indent=2))
