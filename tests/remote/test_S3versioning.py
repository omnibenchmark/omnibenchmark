"""Unit tests for S3versioning module."""

import datetime
from unittest.mock import Mock
import pytest

from omnibenchmark.remote.S3versioning import get_s3_object_versions_and_tags
from omnibenchmark.remote.exception import MinIOStorageBucketManipulationException


@pytest.mark.short
class TestGetS3ObjectVersionsAndTags:
    """Test the get_s3_object_versions_and_tags function."""

    def test_raises_exception_when_bucket_does_not_exist(self):
        """Test that function raises exception when bucket doesn't exist."""
        mock_client = Mock()
        mock_client.bucket_exists.return_value = False

        with pytest.raises(MinIOStorageBucketManipulationException) as exc_info:
            get_s3_object_versions_and_tags(mock_client, "nonexistent_benchmark")

        assert "nonexistent_benchmark does not exist" in str(exc_info.value)
        mock_client.bucket_exists.assert_called_once_with("nonexistent_benchmark")

    def test_returns_empty_dict_when_no_objects(self):
        """Test that function returns empty dict when bucket is empty."""
        mock_client = Mock()
        mock_client.bucket_exists.return_value = True
        mock_client.list_objects.return_value = iter([])

        result = get_s3_object_versions_and_tags(mock_client, "test_benchmark")

        assert result == {}
        mock_client.bucket_exists.assert_called_once_with("test_benchmark")
        mock_client.list_objects.assert_called_once_with(
            "test_benchmark",
            include_version=True,
            include_user_meta=True,
            recursive=True,
        )

    def test_returns_object_metadata_without_tags_in_readonly_mode(self):
        """Test that function returns object metadata without fetching tags in readonly mode."""
        mock_client = Mock()
        mock_client.bucket_exists.return_value = True

        # Create mock object with required attributes
        mock_object = Mock()
        mock_object.object_name = "test_file.txt"
        mock_object.version_id = "v1"
        mock_object.size = 1024
        mock_object.last_modified = datetime.datetime(2025, 1, 1, 12, 0, 0)
        mock_object.is_delete_marker = False
        mock_object.etag = "abc123"

        mock_client.list_objects.return_value = iter([mock_object])

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

        # Should not call get_object_tags in readonly mode
        mock_client.get_object_tags.assert_not_called()

    def test_returns_object_metadata_with_tags_in_readwrite_mode(self):
        """Test that function returns object metadata with tags in read-write mode."""
        mock_client = Mock()
        mock_client.bucket_exists.return_value = True

        mock_object = Mock()
        mock_object.object_name = "test_file.txt"
        mock_object.version_id = "v1"
        mock_object.size = 2048
        mock_object.last_modified = datetime.datetime(2025, 1, 2, 14, 30, 0)
        mock_object.is_delete_marker = False
        mock_object.etag = "def456"

        mock_client.list_objects.return_value = iter([mock_object])
        mock_client.get_object_tags.return_value = {
            "0.1": "valid_version",
            "invalid_tag": "should_be_filtered",
        }

        result = get_s3_object_versions_and_tags(
            mock_client, "test_benchmark", readonly=False
        )

        assert "test_file.txt" in result
        assert "v1" in result["test_file.txt"]
        # Should only include valid version tags
        assert result["test_file.txt"]["v1"]["tags"] == {"0.1": "valid_version"}
        assert result["test_file.txt"]["v1"]["size"] == 2048
        assert result["test_file.txt"]["v1"]["etag"] == "def456"

        mock_client.get_object_tags.assert_called_once_with(
            "test_benchmark", "test_file.txt", version_id="v1"
        )

    def test_handles_multiple_versions_of_same_object(self):
        """Test that function handles multiple versions of the same object."""
        mock_client = Mock()
        mock_client.bucket_exists.return_value = True

        # Create two versions of the same file
        mock_object_v1 = Mock()
        mock_object_v1.object_name = "file.txt"
        mock_object_v1.version_id = "v1"
        mock_object_v1.size = 100
        mock_object_v1.last_modified = datetime.datetime(2025, 1, 1, 10, 0, 0)
        mock_object_v1.is_delete_marker = False
        mock_object_v1.etag = "etag1"

        mock_object_v2 = Mock()
        mock_object_v2.object_name = "file.txt"
        mock_object_v2.version_id = "v2"
        mock_object_v2.size = 200
        mock_object_v2.last_modified = datetime.datetime(2025, 1, 2, 10, 0, 0)
        mock_object_v2.is_delete_marker = False
        mock_object_v2.etag = "etag2"

        mock_client.list_objects.return_value = iter([mock_object_v1, mock_object_v2])
        mock_client.get_object_tags.return_value = {"0.1": "tag"}

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
        """Test that function handles multiple different objects."""
        mock_client = Mock()
        mock_client.bucket_exists.return_value = True

        mock_object1 = Mock()
        mock_object1.object_name = "file1.txt"
        mock_object1.version_id = "v1"
        mock_object1.size = 100
        mock_object1.last_modified = datetime.datetime(2025, 1, 1)
        mock_object1.is_delete_marker = False
        mock_object1.etag = "etag1"

        mock_object2 = Mock()
        mock_object2.object_name = "file2.txt"
        mock_object2.version_id = "v1"
        mock_object2.size = 200
        mock_object2.last_modified = datetime.datetime(2025, 1, 2)
        mock_object2.is_delete_marker = False
        mock_object2.etag = "etag2"

        mock_client.list_objects.return_value = iter([mock_object1, mock_object2])
        mock_client.get_object_tags.return_value = {"0.1": "tag"}

        result = get_s3_object_versions_and_tags(
            mock_client, "test_benchmark", readonly=False
        )

        assert "file1.txt" in result
        assert "file2.txt" in result
        assert result["file1.txt"]["v1"]["size"] == 100
        assert result["file2.txt"]["v1"]["size"] == 200

    def test_skips_tag_retrieval_for_delete_markers(self):
        """Test that function does not retrieve tags for delete markers."""
        mock_client = Mock()
        mock_client.bucket_exists.return_value = True

        mock_delete_marker = Mock()
        mock_delete_marker.object_name = "deleted_file.txt"
        mock_delete_marker.version_id = "v1"
        mock_delete_marker.size = 0
        mock_delete_marker.last_modified = datetime.datetime(2025, 1, 1)
        mock_delete_marker.is_delete_marker = True
        mock_delete_marker.etag = ""

        mock_client.list_objects.return_value = iter([mock_delete_marker])

        result = get_s3_object_versions_and_tags(
            mock_client, "test_benchmark", readonly=False
        )

        assert "deleted_file.txt" in result
        assert result["deleted_file.txt"]["v1"]["tags"] == {}
        assert result["deleted_file.txt"]["v1"]["is_delete_marker"] is True
        # Should not call get_object_tags for delete markers
        mock_client.get_object_tags.assert_not_called()

    def test_filters_out_invalid_version_tags(self):
        """Test that function filters out tags with invalid version keys."""
        mock_client = Mock()
        mock_client.bucket_exists.return_value = True

        mock_object = Mock()
        mock_object.object_name = "file.txt"
        mock_object.version_id = "v1"
        mock_object.size = 100
        mock_object.last_modified = datetime.datetime(2025, 1, 1)
        mock_object.is_delete_marker = False
        mock_object.etag = "etag1"

        mock_client.list_objects.return_value = iter([mock_object])
        mock_client.get_object_tags.return_value = {
            "0.1": "valid",
            "1.2.3": "valid",
            "invalid": "invalid",
            "abc": "invalid",
            "": "invalid",
        }

        result = get_s3_object_versions_and_tags(
            mock_client, "test_benchmark", readonly=False
        )

        # Should only include valid version tags
        assert result["file.txt"]["v1"]["tags"] == {"0.1": "valid", "1.2.3": "valid"}

    def test_handles_none_tags(self):
        """Test that function handles None tags gracefully."""
        mock_client = Mock()
        mock_client.bucket_exists.return_value = True

        mock_object = Mock()
        mock_object.object_name = "file.txt"
        mock_object.version_id = "v1"
        mock_object.size = 100
        mock_object.last_modified = datetime.datetime(2025, 1, 1)
        mock_object.is_delete_marker = False
        mock_object.etag = "etag1"

        mock_client.list_objects.return_value = iter([mock_object])
        mock_client.get_object_tags.return_value = None

        result = get_s3_object_versions_and_tags(
            mock_client, "test_benchmark", readonly=False
        )

        assert result["file.txt"]["v1"]["tags"] == {}

    def test_handles_none_object_name(self):
        """Test that function handles None object name gracefully."""
        mock_client = Mock()
        mock_client.bucket_exists.return_value = True

        mock_object = Mock()
        mock_object.object_name = None
        mock_object.version_id = "v1"
        mock_object.size = 100
        mock_object.last_modified = datetime.datetime(2025, 1, 1)
        mock_object.is_delete_marker = False
        mock_object.etag = "etag1"

        mock_client.list_objects.return_value = iter([mock_object])

        result = get_s3_object_versions_and_tags(
            mock_client, "test_benchmark", readonly=False
        )

        # Should not call get_object_tags when object_name is None
        mock_client.get_object_tags.assert_not_called()
        assert None in result

    def test_handles_complex_scenario_with_mixed_objects(self):
        """Test complex scenario with multiple objects, versions, and delete markers."""
        mock_client = Mock()
        mock_client.bucket_exists.return_value = True

        # Create multiple objects with different states
        objects = []

        # File 1 - two versions
        obj1_v1 = Mock()
        obj1_v1.object_name = "file1.txt"
        obj1_v1.version_id = "v1"
        obj1_v1.size = 100
        obj1_v1.last_modified = datetime.datetime(2025, 1, 1)
        obj1_v1.is_delete_marker = False
        obj1_v1.etag = "etag1"
        objects.append(obj1_v1)

        obj1_v2 = Mock()
        obj1_v2.object_name = "file1.txt"
        obj1_v2.version_id = "v2"
        obj1_v2.size = 150
        obj1_v2.last_modified = datetime.datetime(2025, 1, 2)
        obj1_v2.is_delete_marker = False
        obj1_v2.etag = "etag2"
        objects.append(obj1_v2)

        # File 2 - one version with delete marker
        obj2_v1 = Mock()
        obj2_v1.object_name = "file2.txt"
        obj2_v1.version_id = "v1"
        obj2_v1.size = 200
        obj2_v1.last_modified = datetime.datetime(2025, 1, 1)
        obj2_v1.is_delete_marker = False
        obj2_v1.etag = "etag3"
        objects.append(obj2_v1)

        obj2_v2 = Mock()
        obj2_v2.object_name = "file2.txt"
        obj2_v2.version_id = "v2"
        obj2_v2.size = 0
        obj2_v2.last_modified = datetime.datetime(2025, 1, 3)
        obj2_v2.is_delete_marker = True
        obj2_v2.etag = ""
        objects.append(obj2_v2)

        mock_client.list_objects.return_value = iter(objects)

        # Mock get_object_tags to return different tags for different calls
        def get_tags_side_effect(bucket, obj_name, version_id):
            if obj_name == "file1.txt":
                return {"0.1": "tag1", "invalid": "tag2"}
            elif obj_name == "file2.txt" and version_id == "v1":
                return {"0.2": "tag3"}
            return {}

        mock_client.get_object_tags.side_effect = get_tags_side_effect

        result = get_s3_object_versions_and_tags(
            mock_client, "test_benchmark", readonly=False
        )

        # Verify file1.txt has both versions with tags
        assert "file1.txt" in result
        assert "v1" in result["file1.txt"]
        assert "v2" in result["file1.txt"]
        assert result["file1.txt"]["v1"]["tags"] == {"0.1": "tag1"}
        assert result["file1.txt"]["v2"]["tags"] == {"0.1": "tag1"}
        assert result["file1.txt"]["v1"]["size"] == 100
        assert result["file1.txt"]["v2"]["size"] == 150

        # Verify file2.txt has both versions, including delete marker
        assert "file2.txt" in result
        assert "v1" in result["file2.txt"]
        assert "v2" in result["file2.txt"]
        assert result["file2.txt"]["v1"]["tags"] == {"0.2": "tag3"}
        assert result["file2.txt"]["v2"]["tags"] == {}  # Delete marker has no tags
        assert result["file2.txt"]["v2"]["is_delete_marker"] is True
