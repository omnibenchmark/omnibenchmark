"""
Unit tests for the ConfigAccessor class in omnibenchmark.config module.
"""

import os
import pytest
import tempfile
from pathlib import Path

from omnibenchmark.config import ConfigAccessor, config_dir


@pytest.fixture
def temp_config_file():
    """Create a temporary config file for testing."""
    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as f:
        f.write("""
[test]
key1 = value1
key2 = value2

[another_section]
key3 = value3
        """)
        temp_path = f.name

    yield Path(temp_path)

    # Clean up the temporary file after the test
    os.unlink(temp_path)


@pytest.fixture
def empty_config_file():
    """Create an empty temporary config file for testing."""
    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as f:
        temp_path = f.name

    yield Path(temp_path)

    # Clean up the temporary file after the test
    os.unlink(temp_path)


@pytest.mark.short
def test_config_accessor_get_existing(temp_config_file):
    """Test getting existing values from the config."""
    config = ConfigAccessor(temp_config_file)

    assert config.get("test", "key1") == "value1"
    assert config.get("test", "key2") == "value2"
    assert config.get("another_section", "key3") == "value3"


@pytest.mark.short
def test_config_accessor_get_missing_with_default(temp_config_file):
    """Test getting missing values with defaults."""
    config = ConfigAccessor(temp_config_file)

    # Missing key in existing section
    assert config.get("test", "missing_key", default="default") == "default"

    # Missing section
    assert config.get("missing_section", "key", default="default") == "default"


@pytest.mark.short
def test_config_accessor_get_missing_without_default(temp_config_file):
    """Test getting missing values without defaults."""
    config = ConfigAccessor(temp_config_file)

    # Missing key in existing section
    assert config.get("test", "missing_key") is None

    # Missing section
    assert config.get("missing_section", "key") is None


@pytest.mark.short
def test_config_accessor_set_and_save(empty_config_file):
    """Test setting values and saving to config file."""
    config = ConfigAccessor(empty_config_file)

    # Set some values
    config.set("new_section", "new_key", "new_value")
    config.set("another_section", "key1", "value1")

    # Save to file
    config.save()

    # Create a new config accessor to read from the file
    new_config = ConfigAccessor(empty_config_file)

    # Verify the values were saved
    assert new_config.get("new_section", "new_key") == "new_value"
    assert new_config.get("another_section", "key1") == "value1"


@pytest.mark.short
def test_config_accessor_update_existing(temp_config_file):
    """Test updating existing values."""
    config = ConfigAccessor(temp_config_file)

    # Verify original value
    assert config.get("test", "key1") == "value1"

    # Update the value
    config.set("test", "key1", "updated_value")

    # Verify the value was updated in memory
    assert config.get("test", "key1") == "updated_value"

    # Save and verify persistence
    config.save()
    new_config = ConfigAccessor(temp_config_file)
    assert new_config.get("test", "key1") == "updated_value"


@pytest.mark.short
def test_config_accessor_sections(temp_config_file):
    """Test getting all sections."""
    config = ConfigAccessor(temp_config_file)

    sections = config.sections()
    assert "test" in sections
    assert "another_section" in sections
    assert len(sections) == 2


@pytest.mark.short
def test_config_accessor_options(temp_config_file):
    """Test getting options in a section."""
    config = ConfigAccessor(temp_config_file)

    # Options in existing section
    options = config.options("test")
    assert "key1" in options
    assert "key2" in options
    assert len(options) == 2

    # Options in non-existing section should return empty list
    assert config.options("non_existing") == []


@pytest.mark.short
def test_config_accessor_empty_file(empty_config_file):
    """Test working with an empty config file."""
    config = ConfigAccessor(empty_config_file)

    # Empty file should have no sections
    assert len(config.sections()) == 0

    # Should handle gets gracefully
    assert config.get("section", "key", default="default") == "default"


@pytest.mark.short
def test_default_config_path():
    """Test that ConfigAccessor uses the default path when none is provided."""
    # Create a ConfigAccessor with no path
    config = ConfigAccessor()

    # Check that it's using the expected default path
    assert config.config_path == config_dir / "omnibenchmark.cfg"
