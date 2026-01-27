"""
Unit tests for git module configuration functionality.
"""

import pytest
import tempfile
from pathlib import Path
from unittest.mock import patch

from omnibenchmark.config import get_git_modules_dir, ConfigAccessor


@pytest.mark.short
def test_get_git_modules_dir_default():
    """Test that get_git_modules_dir returns the default .omnibenchmark/git directory."""
    git_dir = get_git_modules_dir()

    # Should default to .omnibenchmark/git
    assert git_dir.name == "git"
    assert git_dir.parent.name == ".omnibenchmark"

    # Directory should be created
    assert git_dir.exists()
    assert git_dir.is_dir()


@pytest.mark.short
def test_get_git_modules_dir_custom_config(tmp_path):
    """Test that get_git_modules_dir respects custom configuration."""
    # Create a custom config file
    custom_git_dir = tmp_path / "custom_git_cache"
    config_file = tmp_path / "test_config.cfg"

    # Create config accessor with custom path
    config = ConfigAccessor(config_file)
    config.set("dirs", "git_modules", str(custom_git_dir))
    config.save()

    # Mock the global config to use our custom one
    with patch("omnibenchmark.config.config", config):
        git_dir = get_git_modules_dir()

        assert git_dir == custom_git_dir
        assert git_dir.exists()
        assert git_dir.is_dir()


@pytest.mark.short
def test_get_git_modules_dir_expanduser():
    """Test that get_git_modules_dir properly expands user home directory."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        config_file = Path(tmp_dir) / "test_config.cfg"

        # Create config with ~ in path
        config = ConfigAccessor(config_file)
        config.set("dirs", "git_modules", "~/.test_omnibenchmark_git")
        config.save()

        # Mock the global config to use our custom one
        with patch("omnibenchmark.config.config", config):
            git_dir = get_git_modules_dir()

            # Should expand ~ to actual home directory
            assert not str(git_dir).startswith("~")
            assert git_dir.name == ".test_omnibenchmark_git"

            # Clean up
            if git_dir.exists():
                git_dir.rmdir()


@pytest.mark.short
def test_config_accessor_git_modules_fallback():
    """Test that ConfigAccessor falls back to default when git_modules key doesn't exist."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        config_file = Path(tmp_dir) / "empty_config.cfg"

        # Create empty config
        config = ConfigAccessor(config_file)

        # Should return default when key doesn't exist
        result = config.get("dirs", "git_modules", ".omnibenchmark/git")
        assert result == ".omnibenchmark/git"


@pytest.mark.short
def test_config_accessor_sections_and_options():
    """Test ConfigAccessor sections and options methods."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        config_file = Path(tmp_dir) / "test_config.cfg"

        config = ConfigAccessor(config_file)
        config.set("dirs", "git_modules", "/custom/git/path")
        config.set("dirs", "datasets", "/custom/datasets/path")
        config.set("other", "setting", "value")
        config.save()

        # Test sections
        sections = config.sections()
        assert "dirs" in sections
        assert "other" in sections

        # Test options in dirs section
        options = config.options("dirs")
        assert "git_modules" in options
        assert "datasets" in options

        # Test options in non-existent section
        empty_options = config.options("nonexistent")
        assert empty_options == []
