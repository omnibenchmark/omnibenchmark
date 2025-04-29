import os
import pytest
import shutil
from pathlib import Path

from omnibenchmark.benchmark.params import Params
from omnibenchmark.benchmark.symlinks import SymlinkManager


@pytest.fixture
def test_dir(tmp_path):
    """Provide a clean temporary directory for each test."""
    yield tmp_path
    # Cleanup after test
    shutil.rmtree(tmp_path)


@pytest.fixture
def manager(test_dir):
    """Provide a SymlinkManager instance for testing."""
    return SymlinkManager(base_dir=test_dir)


@pytest.mark.short
def test_basic_storage(manager):
    """Test basic parameter storage and retrieval."""
    params = Params({"a": 1, "b": 2})
    info = manager.store(params)

    # Check that both directories exist
    assert (Path(manager.base_dir) / info["folder"]).exists()
    # The human name should now be "a_1-b_2"
    assert (Path(manager.base_dir) / info["human"]).exists()

    # Check that we can retrieve the params
    retrieved = manager.get_params(info["human"])
    assert retrieved == params


@pytest.mark.short
def test_human_readable_naming(manager):
    """Test the human-readable naming format."""
    test_cases = [
        ({"a": 1, "b": 2}, "a-1_b-2"),
        ({"foo": "bar"}, "foo-bar"),
        ({"x": 1, "y": 2, "z": 3}, "x-1_y-2_z-3"),
        ({"with space": "has space"}, "with_space-has_space"),
    ]

    for params_dict, expected_name in test_cases:
        params = Params(params_dict)
        info = manager.store(params)
        assert info["human"] == expected_name

        # Verify retrieval works
        retrieved = manager.get_params(info["human"])
        assert retrieved == params


@pytest.mark.short
def test_long_name_fallback(manager):
    """Test that very long parameter names fall back to hash."""
    # Create params with long values
    long_params = Params({"key": "a" * 300})  # This should create a too-long name

    info = manager.store(long_params)
    assert len(info["human"]) <= 255

    # Should still be retrievable
    retrieved = manager.get_params(info["human"])
    assert retrieved == long_params


@pytest.mark.short
def test_invalid_retrieval(manager):
    """Test error handling for non-existent parameters."""
    with pytest.raises(FileNotFoundError):
        manager.get_params("nonexistent_params")


@pytest.mark.short
def test_update_existing(manager):
    """Test updating parameters with same human-readable name."""
    params1 = Params({"a": 1})
    params2 = Params({"a": 2})

    _ = manager.store(params1)
    info2 = manager.store(params2)

    # Should get the updated params
    retrieved = manager.get_params(info2["human"])
    assert retrieved == params2


@pytest.mark.short
def test_multiple_params(manager):
    """Test storing multiple parameter sets."""
    params_list = [Params({"a": 1}), Params({"b": 2}), Params({"c": 3})]

    stored_infos = []
    for params in params_list:
        info = manager.store(params)
        stored_infos.append(info)

    # All should be retrievable
    for info in stored_infos:
        retrieved = manager.get_params(info["human"])
        assert retrieved == info["params"]


@pytest.mark.short
def test_symlink_points_to_correct_location(manager):
    """Test that symlinks point to the correct storage location."""
    params = Params({"a": 1})
    info = manager.store(params)

    symlink_path = Path(manager.base_dir) / info["human"]
    target = symlink_path.readlink()

    # Resolve the target relative to the symlink location
    actual_target = (symlink_path.parent / target).resolve()
    expected_target = Path(manager.base_dir / info["folder"]).resolve()

    assert actual_target == expected_target


@pytest.mark.short
def test_symlink_relative_path(manager):
    """Test that symlinks use correct relative paths."""
    params = Params({"a": 1})
    info = manager.store(params)

    symlink_path = Path(manager.base_dir) / info["human"]

    # Read the actual symlink target (before resolution)
    target = os.readlink(symlink_path)

    # Should point to correct file
    assert (symlink_path.parent / target).resolve() == (
        Path(manager.base_dir) / info["folder"]
    ).resolve()

    # Actually try to read the params file
    params_file = (symlink_path.parent / target) / "parameters.json"
    assert params_file.exists()
    with open(params_file) as f:
        data = f.read()
    assert Params.deserialize(data) == params


@pytest.mark.short
def test_human_folder_identifier_single_param(manager):
    """Test folder identifier with a single parameter."""
    p = Params({"key1": "value1"})
    info = manager.store(p)
    expected_folder_id = "key1-value1"
    assert info["human"] == expected_folder_id


@pytest.mark.short
def test_human_folder_identifier_multiple_params(manager):
    """Test folder identifier with multiple parameters."""
    p = Params({"key1": "value1", "key2": "value2"})
    info = manager.store(p)
    expected_folder_id = "key1-value1_key2-value2"
    assert info["human"] == expected_folder_id


@pytest.mark.short
def test_human_folder_identifier_max_len_exceeded(manager):
    """Test folder identifier with max length exceeded."""
    p = Params({"key1": "value1" * 50, "key2": "value2" * 50})
    info = manager.store(p)
    assert info["human"].endswith("_" + p.hash()[:8])


@pytest.mark.short
def test_human_folder_identifier_with_numeric_values(manager):
    """Test folder identifier with numeric values in keys and values."""
    p = Params({"123": 456, "789": 101112})
    info = manager.store(p)
    expected_folder_id = "123-456_789-101112"
    assert info["human"] == expected_folder_id


@pytest.mark.short
def test_human_folder_identifier_empty_key_value(manager):
    """Test folder identifier with empty key-value pairs."""
    p = Params({"key1": "", "key2": ""})
    info = manager.store(p)
    expected_folder_id = "key1-_key2-"
    assert info["human"] == expected_folder_id


@pytest.mark.short
def test_human_folder_identifier_with_boolean_values(manager):
    """Test folder identifier with boolean values."""
    p = Params({"key1": True, "key2": False})
    info = manager.store(p)
    expected_folder_id = "key1-True_key2-False"
    assert info["human"] == expected_folder_id
