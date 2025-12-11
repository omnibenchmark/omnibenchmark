"""Behavioral unit tests for versioning module."""

import pytest
from datetime import datetime

from omnibenchmark.remote.versioning import (
    get_objects_to_tag,
    filter_objects_to_tag,
    get_single_remoteversion_from_bmversion,
    get_remoteversion_from_bmversion,
    prepare_csv_remoteversion_from_bmversion,
)
from omnibenchmark.remote.exception import MinIOStorageVersioningCorruptionException


class FakeStorageOptions:
    """Minimal storage options for testing."""

    def __init__(self, tracked_dirs=None, results_dirs=None, extra_files=None):
        self.tracked_directories = tracked_dirs or []
        self.results_directories = results_dirs or []
        self.extra_files_to_version_not_in_benchmark_yaml = extra_files or []


class FakeBenchmark:
    """Minimal benchmark for testing."""

    def __init__(self, output_paths=None):
        self._output_paths = output_paths or set()

    def get_output_paths(self):
        return self._output_paths


@pytest.mark.short
class TestGetObjectsToTag:
    """Test get_objects_to_tag function behavior."""

    def test_returns_empty_lists_when_no_objects(self):
        """Should return empty lists when there are no objects."""
        storage_options = FakeStorageOptions(tracked_dirs=["data"])

        result_names, result_versions = get_objects_to_tag({}, storage_options)

        assert result_names == []
        assert result_versions == []

    def test_tags_newest_version_of_object(self):
        """Should select the newest version based on last_modified timestamp."""
        objdic = {
            "data/file.txt": {
                "v1": {
                    "last_modified": datetime(2025, 1, 1),
                    "is_delete_marker": False,
                },
                "v2": {
                    "last_modified": datetime(2025, 1, 2),
                    "is_delete_marker": False,
                },
            }
        }
        storage_options = FakeStorageOptions(tracked_dirs=["data"])

        result_names, result_versions = get_objects_to_tag(objdic, storage_options)

        assert result_names == ["data/file.txt"]
        assert result_versions == ["v2"]

    def test_skips_delete_markers(self):
        """Should skip objects where newest version is a delete marker."""
        objdic = {
            "data/file.txt": {
                "v1": {
                    "last_modified": datetime(2025, 1, 1),
                    "is_delete_marker": False,
                },
                "v2": {
                    "last_modified": datetime(2025, 1, 2),
                    "is_delete_marker": True,
                },
            }
        }
        storage_options = FakeStorageOptions(tracked_dirs=["data"])

        result_names, result_versions = get_objects_to_tag(objdic, storage_options)

        assert result_names == []
        assert result_versions == []

    def test_only_tags_objects_in_tracked_directories(self):
        """Should only tag objects in tracked directories."""
        objdic = {
            "data/file1.txt": {
                "v1": {
                    "last_modified": datetime(2025, 1, 1),
                    "is_delete_marker": False,
                },
            },
            "untracked/file2.txt": {
                "v1": {
                    "last_modified": datetime(2025, 1, 1),
                    "is_delete_marker": False,
                },
            },
        }
        storage_options = FakeStorageOptions(tracked_dirs=["data"])

        result_names, result_versions = get_objects_to_tag(objdic, storage_options)

        assert result_names == ["data/file1.txt"]
        assert result_versions == ["v1"]

    def test_handles_nested_paths(self):
        """Should handle nested directory paths correctly."""
        objdic = {
            "data/subdir/file.txt": {
                "v1": {
                    "last_modified": datetime(2025, 1, 1),
                    "is_delete_marker": False,
                },
            }
        }
        storage_options = FakeStorageOptions(tracked_dirs=["data"])

        result_names, result_versions = get_objects_to_tag(objdic, storage_options)

        assert result_names == ["data/subdir/file.txt"]
        assert result_versions == ["v1"]

    def test_sorts_objects_alphabetically(self):
        """Should return objects in sorted order."""
        objdic = {
            "data/z.txt": {
                "v1": {
                    "last_modified": datetime(2025, 1, 1),
                    "is_delete_marker": False,
                },
            },
            "data/a.txt": {
                "v1": {
                    "last_modified": datetime(2025, 1, 1),
                    "is_delete_marker": False,
                },
            },
            "data/m.txt": {
                "v1": {
                    "last_modified": datetime(2025, 1, 1),
                    "is_delete_marker": False,
                },
            },
        }
        storage_options = FakeStorageOptions(tracked_dirs=["data"])

        result_names, result_versions = get_objects_to_tag(objdic, storage_options)

        assert result_names == ["data/a.txt", "data/m.txt", "data/z.txt"]


@pytest.mark.short
class TestGetExpectedBenchmarkOutputFiles:
    """Test get_expected_benchmark_output_files behavior via filter_objects_to_tag."""

    def test_includes_extra_files_from_glob(self, tmp_path):
        """Should include extra files matched by glob expressions."""
        from omnibenchmark.remote.versioning import get_expected_benchmark_output_files

        # Create some test files
        (tmp_path / "config.yaml").write_text("test")
        (tmp_path / "README.md").write_text("test")

        import os

        original_cwd = os.getcwd()
        os.chdir(tmp_path)

        try:
            storage_options = FakeStorageOptions(extra_files=["*.yaml", "README.md"])
            benchmark = FakeBenchmark(output_paths={"data/output.txt"})

            result = get_expected_benchmark_output_files(benchmark, storage_options)

            # Should include benchmark outputs plus extra files
            assert "data/output.txt" in result
            assert "config.yaml" in result
            assert "README.md" in result
        finally:
            os.chdir(original_cwd)


@pytest.mark.short
class TestFilterObjectsToTag:
    """Test filter_objects_to_tag function behavior."""

    def test_returns_all_objects_when_no_benchmark(self):
        """Should return all objects when benchmark is None."""
        objects = ["data/file1.txt", "data/file2.txt"]
        versions = ["v1", "v2"]
        storage_options = FakeStorageOptions()

        result_objs, result_vers = filter_objects_to_tag(
            objects, versions, storage_options, benchmark=None
        )

        assert result_objs == objects
        assert result_vers == versions

    def test_filters_based_on_benchmark_outputs(self):
        """Should keep only objects that are in benchmark outputs."""
        objects = ["data/output1.txt", "data/output2.txt", "data/other.txt"]
        versions = ["v1", "v2", "v3"]
        storage_options = FakeStorageOptions(
            tracked_dirs=["data"], results_dirs=["data"]
        )
        benchmark = FakeBenchmark(output_paths={"data/output1.txt", "data/output2.txt"})

        result_objs, result_vers = filter_objects_to_tag(
            objects, versions, storage_options, benchmark
        )

        assert set(result_objs) == {"data/output1.txt", "data/output2.txt"}
        assert len(result_vers) == 2

    def test_keeps_objects_in_non_results_directories(self):
        """Should keep all objects in tracked directories that are not results directories."""
        objects = ["config/settings.yaml", "data/output.txt"]
        versions = ["v1", "v2"]
        storage_options = FakeStorageOptions(
            tracked_dirs=["config", "data"], results_dirs=["data"]
        )
        benchmark = FakeBenchmark(output_paths={"data/output.txt"})

        result_objs, result_vers = filter_objects_to_tag(
            objects, versions, storage_options, benchmark
        )

        # config/settings.yaml should be kept even though not in benchmark outputs
        assert "config/settings.yaml" in result_objs
        assert "data/output.txt" in result_objs

    def test_returns_empty_when_no_matches(self):
        """Should return empty lists when no objects match criteria."""
        objects = ["data/unwanted.txt"]
        versions = ["v1"]
        storage_options = FakeStorageOptions(
            tracked_dirs=["data"], results_dirs=["data"]
        )
        benchmark = FakeBenchmark(output_paths={"data/wanted.txt"})

        result_objs, result_vers = filter_objects_to_tag(
            objects, versions, storage_options, benchmark
        )

        assert result_objs == []
        assert result_vers == []


@pytest.mark.short
class TestGetSingleRemoteversionFromBmversion:
    """Test get_single_remoteversion_from_bmversion function behavior."""

    def test_returns_version_with_matching_tag(self):
        """Should return version ID that has the matching tag."""
        objdic = {
            "file.txt": {
                "v1": {"tags": {"0.1": "valid"}},
                "v2": {"tags": {"0.2": "valid"}},
            }
        }

        result = get_single_remoteversion_from_bmversion(objdic, "file.txt", "0.1")

        assert result == "v1"

    def test_returns_none_when_no_matching_tag(self):
        """Should return None when no version has the matching tag."""
        objdic = {
            "file.txt": {
                "v1": {"tags": {"0.1": "valid"}},
                "v2": {"tags": {"0.2": "valid"}},
            }
        }

        result = get_single_remoteversion_from_bmversion(objdic, "file.txt", "0.3")

        assert result is None

    def test_raises_exception_when_multiple_versions_have_same_tag(self):
        """Should raise exception when multiple versions have the same tag (corruption)."""
        objdic = {
            "file.txt": {
                "v1": {"tags": {"0.1": "valid"}},
                "v2": {"tags": {"0.1": "valid"}},
            }
        }

        with pytest.raises(MinIOStorageVersioningCorruptionException) as exc_info:
            get_single_remoteversion_from_bmversion(objdic, "file.txt", "0.1")

        assert "Multiple versions found" in str(exc_info.value)
        assert "file.txt" in str(exc_info.value)


@pytest.mark.short
class TestGetRemoteversionFromBmversion:
    """Test get_remoteversion_from_bmversion function behavior."""

    def test_returns_empty_list_when_no_objects(self):
        """Should return empty list when dictionary is empty."""
        result = get_remoteversion_from_bmversion({}, "0.1")

        assert result == []

    def test_returns_metadata_for_matching_versions(self):
        """Should return complete metadata for versions with matching tags."""
        objdic = {
            "file1.txt": {
                "v1": {
                    "tags": {"0.1": "valid"},
                    "last_modified": datetime(2025, 1, 1),
                    "size": 100,
                    "etag": "abc123",
                }
            },
            "file2.txt": {
                "v1": {
                    "tags": {"0.1": "valid"},
                    "last_modified": datetime(2025, 1, 2),
                    "size": 200,
                    "etag": "def456",
                }
            },
        }

        result = get_remoteversion_from_bmversion(objdic, "0.1")

        assert len(result) == 2
        assert result[0] == ["file1.txt", "v1", datetime(2025, 1, 1), 100, "abc123"]
        assert result[1] == ["file2.txt", "v1", datetime(2025, 1, 2), 200, "def456"]

    def test_skips_objects_without_matching_tag(self):
        """Should skip objects that don't have the matching tag."""
        objdic = {
            "file1.txt": {
                "v1": {
                    "tags": {"0.1": "valid"},
                    "last_modified": datetime(2025, 1, 1),
                    "size": 100,
                    "etag": "abc",
                }
            },
            "file2.txt": {
                "v1": {
                    "tags": {"0.2": "valid"},
                    "last_modified": datetime(2025, 1, 2),
                    "size": 200,
                    "etag": "def",
                }
            },
        }

        result = get_remoteversion_from_bmversion(objdic, "0.1")

        assert len(result) == 1
        assert result[0][0] == "file1.txt"

    def test_handles_multiple_versions_per_object(self):
        """Should handle objects with multiple versions correctly."""
        objdic = {
            "file.txt": {
                "v1": {
                    "tags": {"0.2": "valid"},
                    "last_modified": datetime(2025, 1, 1),
                    "size": 100,
                    "etag": "abc",
                },
                "v2": {
                    "tags": {"0.1": "valid"},
                    "last_modified": datetime(2025, 1, 2),
                    "size": 150,
                    "etag": "def",
                },
            }
        }

        result = get_remoteversion_from_bmversion(objdic, "0.1")

        assert len(result) == 1
        assert result[0] == ["file.txt", "v2", datetime(2025, 1, 2), 150, "def"]


@pytest.mark.short
class TestPrepareCsvRemoteversionFromBmversion:
    """Test prepare_csv_remoteversion_from_bmversion function behavior."""

    def test_returns_header_only_when_empty_list(self):
        """Should return only CSV header when list is empty."""
        result = prepare_csv_remoteversion_from_bmversion([])

        assert result == "name,version_id,last_modified,size,etag\n"

    def test_formats_single_entry_correctly(self):
        """Should format a single entry as CSV correctly."""
        summary = [["file.txt", "v1", datetime(2025, 1, 1), 100, "abc123"]]

        result = prepare_csv_remoteversion_from_bmversion(summary)

        expected = "name,version_id,last_modified,size,etag\n"
        expected += "file.txt,v1,2025-01-01 00:00:00,100,abc123\n"
        assert result == expected

    def test_formats_multiple_entries_correctly(self):
        """Should format multiple entries as CSV correctly."""
        summary = [
            ["file1.txt", "v1", datetime(2025, 1, 1, 10, 30), 100, "abc"],
            ["file2.txt", "v2", datetime(2025, 1, 2, 14, 45), 200, "def"],
        ]

        result = prepare_csv_remoteversion_from_bmversion(summary)

        lines = result.split("\n")
        assert lines[0] == "name,version_id,last_modified,size,etag"
        assert lines[1] == "file1.txt,v1,2025-01-01 10:30:00,100,abc"
        assert lines[2] == "file2.txt,v2,2025-01-02 14:45:00,200,def"
        assert lines[3] == ""  # Trailing newline
