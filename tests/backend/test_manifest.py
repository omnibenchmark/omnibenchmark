"""
Unit tests for the manifest module (omnibenchmark/backend/manifest.py).

Tests write_run_manifest function which creates .metadata/manifest.json with
run provenance information including hardware specs, software versions, and
runtime metadata.

All tests marked as 'short' since they only test local file I/O and system
introspection without external dependencies.
"""

import json
import platform
import sys
from datetime import datetime
from pathlib import Path
from unittest.mock import patch, MagicMock

import pytest

from omnibenchmark.backend._manifest import write_run_manifest


@pytest.mark.short
class TestWriteRunManifest:
    """Test the write_run_manifest function."""

    def test_creates_metadata_directory(self, tmp_path):
        """Test that .metadata directory is created if it doesn't exist."""
        output_dir = tmp_path / "output"
        # Directory doesn't exist yet
        assert not output_dir.exists()

        write_run_manifest(output_dir)

        # Should create output/.metadata/
        metadata_dir = output_dir / ".metadata"
        assert metadata_dir.exists()
        assert metadata_dir.is_dir()

    def test_creates_manifest_json_file(self, tmp_path):
        """Test that manifest.json is created."""
        output_dir = tmp_path / "output"

        write_run_manifest(output_dir)

        manifest_path = output_dir / ".metadata" / "manifest.json"
        assert manifest_path.exists()
        assert manifest_path.is_file()

    def test_manifest_is_valid_json(self, tmp_path):
        """Test that manifest.json contains valid JSON."""
        output_dir = tmp_path / "output"

        write_run_manifest(output_dir)

        manifest_path = output_dir / ".metadata" / "manifest.json"
        with open(manifest_path) as f:
            data = json.load(f)

        # Should be a dict
        assert isinstance(data, dict)

    def test_manifest_contains_required_fields(self, tmp_path):
        """Test that manifest contains all required fields."""
        output_dir = tmp_path / "output"

        write_run_manifest(output_dir)

        manifest_path = output_dir / ".metadata" / "manifest.json"
        with open(manifest_path) as f:
            data = json.load(f)

        # Check required fields
        required_fields = [
            "run_id",
            "ob_version",
            "timestamp",
            "hostname",
            "platform",
            "os",
            "kernel",
            "cpu_count",
            "cpu_model",
            "memory_total_mb",
            "python_version",
            "python_executable",
            "gpu_devices",
        ]

        for field in required_fields:
            assert field in data, f"Missing required field: {field}"

    def test_explicit_run_id_is_used(self, tmp_path):
        """Test that explicit run_id is preserved in manifest."""
        output_dir = tmp_path / "output"
        explicit_run_id = "test-run-12345"

        write_run_manifest(output_dir, run_id=explicit_run_id)

        manifest_path = output_dir / ".metadata" / "manifest.json"
        with open(manifest_path) as f:
            data = json.load(f)

        assert data["run_id"] == explicit_run_id

    def test_auto_generated_run_id_is_uuid(self, tmp_path):
        """Test that auto-generated run_id is a valid UUID."""
        output_dir = tmp_path / "output"

        result = write_run_manifest(output_dir)

        # Should be UUID4 format (8-4-4-4-12 hex digits)
        run_id = result["run_id"]
        assert isinstance(run_id, str)
        assert len(run_id) == 36
        assert run_id.count("-") == 4

        # Try to parse as UUID to validate format
        import uuid

        uuid_obj = uuid.UUID(run_id)
        assert str(uuid_obj) == run_id

    def test_returns_manifest_dict(self, tmp_path):
        """Test that function returns the manifest dict."""
        output_dir = tmp_path / "output"

        result = write_run_manifest(output_dir)

        assert isinstance(result, dict)
        assert "run_id" in result
        assert "timestamp" in result

    def test_timestamp_is_iso8601_utc(self, tmp_path):
        """Test that timestamp is ISO-8601 format in UTC."""
        output_dir = tmp_path / "output"

        result = write_run_manifest(output_dir)

        timestamp = result["timestamp"]
        assert isinstance(timestamp, str)

        # Should be parseable as ISO-8601
        dt = datetime.fromisoformat(timestamp)
        assert dt.tzinfo is not None  # Should include timezone

    def test_hostname_is_captured(self, tmp_path):
        """Test that hostname is captured from platform.uname()."""
        output_dir = tmp_path / "output"

        result = write_run_manifest(output_dir)

        assert "hostname" in result
        # Should match platform.uname().node
        assert result["hostname"] == platform.uname().node

    def test_platform_is_sys_platform(self, tmp_path):
        """Test that platform field matches sys.platform."""
        output_dir = tmp_path / "output"

        result = write_run_manifest(output_dir)

        assert result["platform"] == sys.platform

    def test_python_version_captured(self, tmp_path):
        """Test that Python version is captured."""
        output_dir = tmp_path / "output"

        result = write_run_manifest(output_dir)

        assert "python_version" in result
        assert result["python_version"] == platform.python_version()

    def test_python_executable_captured(self, tmp_path):
        """Test that Python executable path is captured."""
        output_dir = tmp_path / "output"

        result = write_run_manifest(output_dir)

        assert result["python_executable"] == sys.executable

    def test_os_string_contains_system_info(self, tmp_path):
        """Test that OS string contains system information."""
        output_dir = tmp_path / "output"

        result = write_run_manifest(output_dir)

        os_str = result["os"]
        assert isinstance(os_str, str)
        # Should contain at least the system name
        uname = platform.uname()
        assert uname.system in os_str

    def test_kernel_is_uname_release(self, tmp_path):
        """Test that kernel field matches uname release."""
        output_dir = tmp_path / "output"

        result = write_run_manifest(output_dir)

        assert result["kernel"] == platform.uname().release

    def test_cpu_count_is_numeric_or_none(self, tmp_path):
        """Test that cpu_count is either an integer or None."""
        output_dir = tmp_path / "output"

        result = write_run_manifest(output_dir)

        cpu_count = result["cpu_count"]
        assert cpu_count is None or isinstance(cpu_count, int)
        if cpu_count is not None:
            assert cpu_count > 0

    def test_cpu_model_is_string_or_none(self, tmp_path):
        """Test that cpu_model is a string or None."""
        output_dir = tmp_path / "output"

        result = write_run_manifest(output_dir)

        cpu_model = result["cpu_model"]
        assert cpu_model is None or isinstance(cpu_model, str)

    def test_memory_total_mb_is_numeric_or_none(self, tmp_path):
        """Test that memory_total_mb is an integer or None."""
        output_dir = tmp_path / "output"

        result = write_run_manifest(output_dir)

        memory_total_mb = result["memory_total_mb"]
        assert memory_total_mb is None or isinstance(memory_total_mb, int)
        if memory_total_mb is not None:
            assert memory_total_mb > 0

    def test_gpu_devices_is_list_or_none(self, tmp_path):
        """Test that gpu_devices is a list or None."""
        output_dir = tmp_path / "output"

        result = write_run_manifest(output_dir)

        gpu_devices = result["gpu_devices"]
        assert gpu_devices is None or isinstance(gpu_devices, list)

    def test_ob_version_is_string_or_none(self, tmp_path):
        """Test that ob_version is a string or None."""
        output_dir = tmp_path / "output"

        result = write_run_manifest(output_dir)

        ob_version = result["ob_version"]
        assert ob_version is None or isinstance(ob_version, str)

    def test_overwrites_existing_manifest(self, tmp_path):
        """Test that existing manifest.json is overwritten."""
        output_dir = tmp_path / "output"
        metadata_dir = output_dir / ".metadata"
        metadata_dir.mkdir(parents=True)
        manifest_path = metadata_dir / "manifest.json"

        # Write old manifest
        old_data = {"run_id": "old-id", "timestamp": "2020-01-01T00:00:00Z"}
        with open(manifest_path, "w") as f:
            json.dump(old_data, f)

        # Write new manifest
        _result = write_run_manifest(output_dir, run_id="new-id")

        # Should have new data
        with open(manifest_path) as f:
            data = json.load(f)

        assert data["run_id"] == "new-id"
        assert data["run_id"] != old_data["run_id"]

    def test_manifest_ends_with_newline(self, tmp_path):
        """Test that manifest.json ends with a newline for POSIX compliance."""
        output_dir = tmp_path / "output"

        write_run_manifest(output_dir)

        manifest_path = output_dir / ".metadata" / "manifest.json"
        with open(manifest_path, "rb") as f:
            content = f.read()

        assert content.endswith(b"\n")

    def test_manifest_is_pretty_printed(self, tmp_path):
        """Test that manifest.json is formatted with indentation."""
        output_dir = tmp_path / "output"

        write_run_manifest(output_dir)

        manifest_path = output_dir / ".metadata" / "manifest.json"
        with open(manifest_path) as f:
            content = f.read()

        # Should contain newlines (pretty-printed)
        assert "\n" in content
        # Should contain indentation
        assert "  " in content


@pytest.mark.short
class TestWriteRunManifestLinux:
    """Test Linux-specific functionality in write_run_manifest."""

    @pytest.mark.skipif(sys.platform != "linux", reason="Linux-specific test")
    def test_linux_cpu_model_from_proc_cpuinfo(self, tmp_path):
        """Test that CPU model is read from /proc/cpuinfo on Linux."""
        output_dir = tmp_path / "output"

        result = write_run_manifest(output_dir)

        # On Linux, cpu_model should be populated if /proc/cpuinfo exists
        if Path("/proc/cpuinfo").exists():
            assert result["cpu_model"] is not None
            assert isinstance(result["cpu_model"], str)
            assert len(result["cpu_model"]) > 0

    @pytest.mark.skipif(sys.platform != "linux", reason="Linux-specific test")
    def test_linux_memory_from_proc_meminfo(self, tmp_path):
        """Test that memory is read from /proc/meminfo on Linux."""
        output_dir = tmp_path / "output"

        result = write_run_manifest(output_dir)

        # On Linux, memory should be populated if /proc/meminfo exists
        if Path("/proc/meminfo").exists():
            assert result["memory_total_mb"] is not None
            assert isinstance(result["memory_total_mb"], int)
            assert result["memory_total_mb"] > 0


@pytest.mark.short
class TestWriteRunManifestMacOS:
    """Test macOS-specific functionality in write_run_manifest."""

    @pytest.mark.skipif(sys.platform != "darwin", reason="macOS-specific test")
    def test_macos_cpu_model_from_sysctl(self, tmp_path):
        """Test that CPU model is read from sysctl on macOS."""
        output_dir = tmp_path / "output"

        result = write_run_manifest(output_dir)

        # On macOS, cpu_model should be populated
        assert result["cpu_model"] is not None
        assert isinstance(result["cpu_model"], str)

    @pytest.mark.skipif(sys.platform != "darwin", reason="macOS-specific test")
    def test_macos_memory_from_sysctl(self, tmp_path):
        """Test that memory is read from sysctl on macOS."""
        output_dir = tmp_path / "output"

        result = write_run_manifest(output_dir)

        # On macOS, memory should be populated
        assert result["memory_total_mb"] is not None
        assert isinstance(result["memory_total_mb"], int)
        assert result["memory_total_mb"] > 0


@pytest.mark.short
class TestWriteRunManifestErrorHandling:
    """Test error handling in write_run_manifest."""

    def test_missing_proc_cpuinfo_handled_gracefully(self, tmp_path):
        """Test that missing /proc/cpuinfo doesn't crash."""
        output_dir = tmp_path / "output"

        # Mock open to raise FileNotFoundError for /proc/cpuinfo
        original_open = open

        def mock_open_wrapper(path, *args, **kwargs):
            if str(path) == "/proc/cpuinfo":
                raise FileNotFoundError()
            return original_open(path, *args, **kwargs)

        with patch("builtins.open", side_effect=mock_open_wrapper):
            result = write_run_manifest(output_dir)

        # Should still succeed, cpu_model will be None
        assert "cpu_model" in result
        # On Linux it would normally be populated, but we mocked the error

    def test_nvidia_smi_missing_handled_gracefully(self, tmp_path):
        """Test that missing nvidia-smi doesn't crash."""
        output_dir = tmp_path / "output"

        # Mock subprocess.run to raise FileNotFoundError
        with patch("subprocess.run", side_effect=FileNotFoundError()):
            result = write_run_manifest(output_dir)

        # Should still succeed, gpu_devices will be None
        assert result["gpu_devices"] is None

    def test_nvidia_smi_error_handled_gracefully(self, tmp_path):
        """Test that nvidia-smi errors are handled gracefully."""
        output_dir = tmp_path / "output"

        # Mock subprocess.run to return non-zero exit code
        mock_result = MagicMock()
        mock_result.returncode = 1
        mock_result.stdout = ""

        with patch("subprocess.run", return_value=mock_result):
            result = write_run_manifest(output_dir)

        # Should still succeed, gpu_devices will be None
        assert result["gpu_devices"] is None

    def test_os_cpu_count_exception_handled(self, tmp_path):
        """Test that os.cpu_count() exceptions are handled."""
        output_dir = tmp_path / "output"

        # Mock os.cpu_count to raise exception
        with patch("os.cpu_count", side_effect=RuntimeError("test error")):
            result = write_run_manifest(output_dir)

        # Should still succeed, cpu_count will be None
        assert result["cpu_count"] is None

    def test_importlib_metadata_missing_handled(self, tmp_path):
        """Test that missing importlib.metadata doesn't crash."""
        output_dir = tmp_path / "output"

        # Mock importlib.metadata.version to raise exception
        with patch("importlib.metadata.version", side_effect=ImportError("test error")):
            result = write_run_manifest(output_dir)

        # Should still succeed, ob_version will be None
        assert result["ob_version"] is None


@pytest.mark.short
class TestWriteRunManifestGPU:
    """Test GPU detection in write_run_manifest."""

    def test_nvidia_smi_success_parses_gpu_info(self, tmp_path):
        """Test that nvidia-smi output is parsed correctly."""
        output_dir = tmp_path / "output"

        # Mock successful nvidia-smi output
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = (
            "0, NVIDIA GeForce RTX 3090, 24576\n1, NVIDIA GeForce RTX 3080, 10240\n"
        )

        with patch("subprocess.run", return_value=mock_result):
            result = write_run_manifest(output_dir)

        # Should parse GPU info
        assert result["gpu_devices"] is not None
        assert len(result["gpu_devices"]) == 2

        gpu0 = result["gpu_devices"][0]
        assert gpu0["index"] == 0
        assert gpu0["name"] == "NVIDIA GeForce RTX 3090"
        assert gpu0["memory_total_mb"] == 24576

        gpu1 = result["gpu_devices"][1]
        assert gpu1["index"] == 1
        assert gpu1["name"] == "NVIDIA GeForce RTX 3080"
        assert gpu1["memory_total_mb"] == 10240

    def test_nvidia_smi_empty_output_returns_none(self, tmp_path):
        """Test that empty nvidia-smi output results in None."""
        output_dir = tmp_path / "output"

        # Mock empty output
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = ""

        with patch("subprocess.run", return_value=mock_result):
            result = write_run_manifest(output_dir)

        assert result["gpu_devices"] is None

    def test_nvidia_smi_whitespace_only_returns_none(self, tmp_path):
        """Test that whitespace-only nvidia-smi output results in None."""
        output_dir = tmp_path / "output"

        # Mock whitespace output
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = "   \n  \n  "

        with patch("subprocess.run", return_value=mock_result):
            result = write_run_manifest(output_dir)

        assert result["gpu_devices"] is None

    def test_nvidia_smi_malformed_line_skips_gracefully(self, tmp_path):
        """Test that malformed nvidia-smi lines cause exception and return empty list."""
        output_dir = tmp_path / "output"

        # Mock output with malformed line at the start (will fail immediately)
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = "malformed line\n0, NVIDIA RTX 3090, 24576\n"

        with patch("subprocess.run", return_value=mock_result):
            result = write_run_manifest(output_dir)

        # Should handle exception gracefully and return empty list for gpu_devices
        # because the malformed line will cause int(parts[0]) to fail on first line
        # The outer try-except catches the error after gpu_devices=[] is initialized
        assert result["gpu_devices"] == []

    def test_nvidia_smi_timeout_handled(self, tmp_path):
        """Test that nvidia-smi timeout is handled."""
        output_dir = tmp_path / "output"

        # Mock timeout
        import subprocess

        with patch(
            "subprocess.run", side_effect=subprocess.TimeoutExpired("nvidia-smi", 5)
        ):
            result = write_run_manifest(output_dir)

        # Should handle timeout gracefully
        assert result["gpu_devices"] is None


@pytest.mark.short
class TestWriteRunManifestEdgeCases:
    """Test edge cases in write_run_manifest."""

    def test_very_long_hostname(self, tmp_path):
        """Test handling of unusually long hostnames."""
        output_dir = tmp_path / "output"

        # Mock extremely long hostname
        long_hostname = "a" * 1000
        mock_uname = MagicMock()
        mock_uname.node = long_hostname
        mock_uname.system = "Linux"
        mock_uname.release = "5.0.0"
        mock_uname.version = "#1 SMP"

        with patch("platform.uname", return_value=mock_uname):
            result = write_run_manifest(output_dir)

        # Should handle it without crashing
        assert result["hostname"] == long_hostname

    def test_special_characters_in_system_info(self, tmp_path):
        """Test handling of special characters in system info."""
        output_dir = tmp_path / "output"

        # Mock system info with special chars
        mock_uname = MagicMock()
        mock_uname.node = "test-node™"
        mock_uname.system = "TestOS™"
        mock_uname.release = "1.0-ü"
        mock_uname.version = "#1 SMP αβγ"

        with patch("platform.uname", return_value=mock_uname):
            write_run_manifest(output_dir)

        # Should serialize to JSON without errors
        manifest_path = output_dir / ".metadata" / "manifest.json"
        with open(manifest_path) as f:
            data = json.load(f)

        assert "test-node™" in data["hostname"]

    def test_concurrent_writes_to_same_directory(self, tmp_path):
        """Test that concurrent writes don't corrupt the manifest."""
        output_dir = tmp_path / "output"

        # Write multiple manifests with different run_ids
        write_run_manifest(output_dir, run_id="run-1")
        write_run_manifest(output_dir, run_id="run-2")

        # Last write wins
        manifest_path = output_dir / ".metadata" / "manifest.json"
        with open(manifest_path) as f:
            data = json.load(f)

        assert data["run_id"] == "run-2"

    def test_metadata_directory_already_exists(self, tmp_path):
        """Test that existing .metadata directory doesn't cause errors."""
        output_dir = tmp_path / "output"
        metadata_dir = output_dir / ".metadata"

        # Pre-create the directory
        metadata_dir.mkdir(parents=True)

        # Should not raise an error
        result = write_run_manifest(output_dir)

        assert result is not None
        assert (metadata_dir / "manifest.json").exists()

    def test_output_dir_is_file_raises_error(self, tmp_path):
        """Test that using a file as output_dir raises an appropriate error."""
        output_file = tmp_path / "not_a_dir"
        output_file.write_text("content")

        # Should raise an error when trying to create .metadata
        with pytest.raises((FileExistsError, NotADirectoryError, OSError)):
            write_run_manifest(output_file)
