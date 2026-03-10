"""Unit tests for S3versioning module."""

import datetime
from unittest.mock import Mock

import pytest
from botocore.exceptions import ClientError

from omnibenchmark.remote.S3versioning import get_s3_object_versions_and_tags
from omnibenchmark.remote.exception import S3StorageBucketManipulationException


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_paginator(*pages):
    """Return a mock boto3 paginator that yields the given pages."""
    paginator = Mock()
    paginator.paginate.return_value = list(pages)
    return paginator


def _version(key, vid, size=100, etag="etag1", last_modified=None):
    """Build a boto3-style Version dict (ETag is quoted as in real responses)."""
    return {
        "Key": key,
        "VersionId": vid,
        "Size": size,
        "ETag": f'"{etag}"',
        "LastModified": last_modified or datetime.datetime(2025, 1, 1),
    }


def _delete_marker(key, vid, last_modified=None):
    """Build a boto3-style DeleteMarker dict."""
    return {
        "Key": key,
        "VersionId": vid,
        "LastModified": last_modified or datetime.datetime(2025, 1, 1),
    }


def _client_error(code="404"):
    return ClientError({"Error": {"Code": code, "Message": "Error"}}, "Op")


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


@pytest.mark.short
class TestGetS3ObjectVersionsAndTags:
    """Test the get_s3_object_versions_and_tags function."""

    def test_raises_exception_when_bucket_does_not_exist(self):
        mock_client = Mock()
        mock_client.head_bucket.side_effect = _client_error("404")

        with pytest.raises(S3StorageBucketManipulationException) as exc_info:
            get_s3_object_versions_and_tags(mock_client, "nonexistent_benchmark")

        assert "nonexistent_benchmark does not exist" in str(exc_info.value)

    def test_returns_empty_dict_when_no_objects(self):
        mock_client = Mock()
        mock_client.get_paginator.return_value = _make_paginator({})

        result = get_s3_object_versions_and_tags(mock_client, "test_benchmark")

        assert result == {}
        mock_client.get_paginator.assert_called_once_with("list_object_versions")
        mock_client.get_paginator.return_value.paginate.assert_called_once_with(
            Bucket="test_benchmark"
        )

    def test_returns_object_metadata_without_tags_in_readonly_mode(self):
        mock_client = Mock()
        mock_client.get_paginator.return_value = _make_paginator(
            {
                "Versions": [
                    _version(
                        "test_file.txt",
                        "v1",
                        size=1024,
                        etag="abc123",
                        last_modified=datetime.datetime(2025, 1, 1, 12, 0, 0),
                    )
                ]
            }
        )

        result = get_s3_object_versions_and_tags(
            mock_client, "test_benchmark", readonly=True
        )

        assert "test_file.txt" in result
        assert "v1" in result["test_file.txt"]
        assert result["test_file.txt"]["v1"]["tags"] == {}
        assert result["test_file.txt"]["v1"]["size"] == 1024
        assert result["test_file.txt"]["v1"]["last_modified"] == datetime.datetime(
            2025, 1, 1, 12, 0, 0
        )
        assert result["test_file.txt"]["v1"]["is_delete_marker"] is False
        assert result["test_file.txt"]["v1"]["etag"] == "abc123"
        # No tagging calls in readonly mode
        mock_client.get_object_tagging.assert_not_called()

    def test_returns_object_metadata_with_tags_in_readwrite_mode(self):
        mock_client = Mock()
        mock_client.get_paginator.return_value = _make_paginator(
            {
                "Versions": [
                    _version(
                        "test_file.txt",
                        "v1",
                        size=2048,
                        etag="def456",
                        last_modified=datetime.datetime(2025, 1, 2, 14, 30, 0),
                    )
                ]
            }
        )
        mock_client.get_object_tagging.return_value = {
            "TagSet": [
                {"Key": "0.1", "Value": "valid_version"},
                {"Key": "invalid_tag", "Value": "should_be_filtered"},
            ]
        }

        result = get_s3_object_versions_and_tags(
            mock_client, "test_benchmark", readonly=False
        )

        assert "test_file.txt" in result
        # Only valid version tags are kept
        assert result["test_file.txt"]["v1"]["tags"] == {"0.1": "valid_version"}
        assert result["test_file.txt"]["v1"]["size"] == 2048
        assert result["test_file.txt"]["v1"]["etag"] == "def456"
        mock_client.get_object_tagging.assert_called_once_with(
            Bucket="test_benchmark", Key="test_file.txt", VersionId="v1"
        )

    def test_handles_multiple_versions_of_same_object(self):
        mock_client = Mock()
        mock_client.get_paginator.return_value = _make_paginator(
            {
                "Versions": [
                    _version("file.txt", "v1", size=100, etag="etag1"),
                    _version("file.txt", "v2", size=200, etag="etag2"),
                ]
            }
        )
        mock_client.get_object_tagging.return_value = {
            "TagSet": [{"Key": "0.1", "Value": "tag"}]
        }

        result = get_s3_object_versions_and_tags(
            mock_client, "test_benchmark", readonly=False
        )

        assert "file.txt" in result
        assert "v1" in result["file.txt"]
        assert "v2" in result["file.txt"]
        assert result["file.txt"]["v1"]["size"] == 100
        assert result["file.txt"]["v2"]["size"] == 200
        assert result["file.txt"]["v1"]["etag"] == "etag1"
        assert result["file.txt"]["v2"]["etag"] == "etag2"

    def test_handles_multiple_different_objects(self):
        mock_client = Mock()
        mock_client.get_paginator.return_value = _make_paginator(
            {
                "Versions": [
                    _version("file1.txt", "v1", size=100, etag="etag1"),
                    _version("file2.txt", "v1", size=200, etag="etag2"),
                ]
            }
        )
        mock_client.get_object_tagging.return_value = {
            "TagSet": [{"Key": "0.1", "Value": "tag"}]
        }

        result = get_s3_object_versions_and_tags(
            mock_client, "test_benchmark", readonly=False
        )

        assert "file1.txt" in result
        assert "file2.txt" in result
        assert result["file1.txt"]["v1"]["size"] == 100
        assert result["file2.txt"]["v1"]["size"] == 200

    def test_skips_tag_retrieval_for_delete_markers(self):
        mock_client = Mock()
        mock_client.get_paginator.return_value = _make_paginator(
            {"DeleteMarkers": [_delete_marker("deleted_file.txt", "v1")]}
        )

        result = get_s3_object_versions_and_tags(
            mock_client, "test_benchmark", readonly=False
        )

        assert "deleted_file.txt" in result
        assert result["deleted_file.txt"]["v1"]["tags"] == {}
        assert result["deleted_file.txt"]["v1"]["is_delete_marker"] is True
        # Delete markers never trigger a tagging call
        mock_client.get_object_tagging.assert_not_called()

    def test_filters_out_invalid_version_tags(self):
        mock_client = Mock()
        mock_client.get_paginator.return_value = _make_paginator(
            {"Versions": [_version("file.txt", "v1")]}
        )
        mock_client.get_object_tagging.return_value = {
            "TagSet": [
                {"Key": "0.1", "Value": "valid"},
                {"Key": "1.2.3", "Value": "valid"},
                {"Key": "invalid", "Value": "invalid"},
                {"Key": "abc", "Value": "invalid"},
                {"Key": "", "Value": "invalid"},
            ]
        }

        result = get_s3_object_versions_and_tags(
            mock_client, "test_benchmark", readonly=False
        )

        assert result["file.txt"]["v1"]["tags"] == {"0.1": "valid", "1.2.3": "valid"}

    def test_silently_skips_tags_on_client_error(self):
        """ClientError from get_object_tagging is silently ignored (e.g. AccessDenied)."""
        mock_client = Mock()
        mock_client.get_paginator.return_value = _make_paginator(
            {"Versions": [_version("file.txt", "v1")]}
        )
        mock_client.get_object_tagging.side_effect = _client_error("AccessDenied")

        result = get_s3_object_versions_and_tags(
            mock_client, "test_benchmark", readonly=False
        )

        assert result["file.txt"]["v1"]["tags"] == {}

    def test_handles_complex_scenario_with_mixed_objects(self):
        """Multiple objects, multiple versions, and delete markers across one page."""
        mock_client = Mock()
        mock_client.get_paginator.return_value = _make_paginator(
            {
                "Versions": [
                    _version(
                        "file1.txt",
                        "v1",
                        size=100,
                        etag="etag1",
                        last_modified=datetime.datetime(2025, 1, 1),
                    ),
                    _version(
                        "file1.txt",
                        "v2",
                        size=150,
                        etag="etag2",
                        last_modified=datetime.datetime(2025, 1, 2),
                    ),
                    _version(
                        "file2.txt",
                        "v1",
                        size=200,
                        etag="etag3",
                        last_modified=datetime.datetime(2025, 1, 1),
                    ),
                ],
                "DeleteMarkers": [
                    _delete_marker(
                        "file2.txt", "v2", last_modified=datetime.datetime(2025, 1, 3)
                    )
                ],
            }
        )

        def _tags(Bucket, Key, VersionId):
            if Key == "file1.txt":
                return {
                    "TagSet": [
                        {"Key": "0.1", "Value": "tag1"},
                        {"Key": "invalid", "Value": "x"},
                    ]
                }
            if Key == "file2.txt" and VersionId == "v1":
                return {"TagSet": [{"Key": "0.2", "Value": "tag3"}]}
            return {"TagSet": []}

        mock_client.get_object_tagging.side_effect = _tags

        result = get_s3_object_versions_and_tags(
            mock_client, "test_benchmark", readonly=False
        )

        assert result["file1.txt"]["v1"]["tags"] == {"0.1": "tag1"}
        assert result["file1.txt"]["v2"]["tags"] == {"0.1": "tag1"}
        assert result["file1.txt"]["v1"]["size"] == 100
        assert result["file1.txt"]["v2"]["size"] == 150

        assert result["file2.txt"]["v1"]["tags"] == {"0.2": "tag3"}
        assert result["file2.txt"]["v2"]["tags"] == {}
        assert result["file2.txt"]["v2"]["is_delete_marker"] is True
