def tests_bucket_policy():
    return {
        "Version": "2012-10-17",
        "Statement": [
            {
                "Effect": "Allow",
                "Action": ["s3:*"],
                "Resource": ["arn:aws:s3:::test.*", "arn:aws:s3:::test?.*"],
            },
            {
                "Effect": "Allow",
                "Action": ["s3:ListAllMyBuckets"],
                "Resource": ["arn:aws:s3:::*"],
            },
        ],
    }


def benchmarker_access_token_policy(benchmark):
    """S3 policy for access token for specific benchmark, allows archiving (Governance Retention)"""
    return {
        "Version": "2012-10-17",
        "Statement": [
            {
                "Effect": "Allow",
                "Action": ["s3:*"],
                "Resource": [
                    f"arn:aws:s3:::{benchmark}.*.*",
                    f"arn:aws:s3:::{benchmark}.*.*/*",
                    f"arn:aws:s3:::{benchmark}.overview/*",
                    f"arn:aws:s3:::{benchmark}.test.?/*",
                ],
            },
            {
                "Effect": "Allow",
                "Action": ["s3:ListAllMyBuckets"],
                "Resource": ["arn:aws:s3:::*"],
            },
            {
                "Effect": "Deny",
                "Action": ["s3:BypassGovernanceRetention"],
                "Resource": [
                    f"arn:aws:s3:::{benchmark}.*.*/*",
                    f"arn:aws:s3:::{benchmark}.overview/*",
                    f"arn:aws:s3:::{benchmark}.test.?/*",
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
                "Action": ["s3:GetBucketLocation", "s3:ListBucket", "s3:ListObjects"],
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


# access key: ekjt8qvAmOeTLRsWcqVL
# secrete key: H0nxqGPB58ZikZuNTpGZPQ9YOtOhofZkvJ0WgDRl
