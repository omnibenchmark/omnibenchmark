import pytest
from omni.benchmark.params import Params


def test_init_empty():
    """Test creating empty Params object."""
    p = Params()
    assert len(p._params) == 0


def test_init_with_dict():
    """Test creating Params from dict."""
    d = {"b": 2, "a": 1}
    p = Params(d)
    assert p["a"] == 1
    assert p["b"] == 2
    # Check if keys are sorted
    assert list(p._params.keys()) == ["a", "b"]


def test_init_with_params():
    """Test creating Params from another Params object."""
    p1 = Params({"b": 2, "a": 1})
    p2 = Params(p1)
    assert p1 == p2


def test_equality():
    """Test equality comparisons."""
    p1 = Params({"b": 2, "a": 1})
    p2 = Params({"a": 1, "b": 2})
    d = {"a": 1, "b": 2}

    assert p1 == p2
    assert p1 == d
    assert p2 == d


def test_hash():
    """Test hash stability."""
    p1 = Params({"b": 2, "a": 1})
    p2 = Params({"a": 1, "b": 2})

    assert p1.hash() == p2.hash()


def test_update():
    """Test updating params."""
    p1 = Params({"a": 1})
    p2 = Params({"b": 2})
    d = {"c": 3}

    p1.update(p2)
    assert p1 == Params({"a": 1, "b": 2})

    p1.update(d)
    assert p1 == Params({"a": 1, "b": 2, "c": 3})


def test_item_access():
    """Test dict-like item access."""
    p = Params({"a": 1})

    assert p["a"] == 1

    p["b"] = 2
    assert p["b"] == 2

    with pytest.raises(KeyError):
        _ = p["nonexistent"]


def test_serialization():
    """Test serialization and deserialization."""
    original = Params({"b": 2, "a": 1})
    serialized = original.serialize()
    deserialized = Params.deserialize(serialized)

    assert original == deserialized
    assert serialized == '{"a": 1, "b": 2}'


def test_to_cli_args_gnu():
    """Test conversion to GNU-style CLI args."""
    p = Params({"foo": "bar", "num": 42})
    args = p.to_cli_args(style="gnu")
    assert args == ["--foo", "bar", "--num", "42"]


def test_to_cli_args_equals():
    """Test conversion to equals-style CLI args."""
    p = Params({"foo": "bar", "num": 42})
    args = p.to_cli_args(style="equals")
    assert args == ["--foo=bar", "--num=42"]


def test_from_cli_args_gnu():
    """Test parsing GNU-style CLI args."""
    args = ["--foo", "bar", "--num", "42"]
    p = Params.from_cli_args(args)
    assert p == Params({"foo": "bar", "num": "42"})


def test_from_cli_args_equals():
    """Test parsing equals-style CLI args."""
    args = ["--foo=bar", "--num=42"]
    p = Params.from_cli_args(args)
    assert p == Params({"foo": "bar", "num": "42"})


def test_from_cli_args_mixed():
    """Test parsing mixed-style CLI args."""
    args = ["--foo=bar", "--num", "42"]
    p = Params.from_cli_args(args)
    assert p == Params({"foo": "bar", "num": "42"})


def test_from_cli_args_with_flags():
    """Test parsing CLI args with flag-like parameters."""
    args = ["--foo", "bar", "--flag", "--num", "42"]
    p = Params.from_cli_args(args)
    assert p == Params({"foo": "bar", "flag": True, "num": "42"})


def test_generate_folder_identifier_empty_params():
    """Test folder identifier with an empty Params object."""
    p = Params()
    folder_id = p.generate_folder_identifier()
    assert folder_id == "", f"Expected empty folder name but got: {folder_id}"


def test_generate_folder_identifier_single_param():
    """Test folder identifier with a single parameter."""
    p = Params({"key1": "value1"})
    folder_id = p.generate_folder_identifier()
    expected_folder_id = "key1-value1"
    assert (
        folder_id == expected_folder_id
    ), f"Expected '{expected_folder_id}' but got: {folder_id}"


def test_generate_folder_identifier_multiple_params():
    """Test folder identifier with multiple parameters."""
    p = Params({"key1": "value1", "key2": "value2"})
    folder_id = p.generate_folder_identifier()
    expected_folder_id = "key1-value1_key2-value2"
    assert (
        folder_id == expected_folder_id
    ), f"Expected '{expected_folder_id}' but got: {folder_id}"


def test_generate_folder_identifier_special_characters():
    """Test folder identifier with special characters in keys and values."""
    p = Params({"key@1": "value#1", "key_2": "value&2"})
    folder_id = p.generate_folder_identifier()
    expected_folder_id = "key1-value1_key_2-value2"
    assert (
        folder_id == expected_folder_id
    ), f"Expected '{expected_folder_id}' but got: {folder_id}"


def test_generate_folder_identifier_max_len_exact():
    """Test folder identifier with max length exactly fitting."""
    p = Params({"key1": "value1" * 10, "key2": "value2" * 10})
    folder_id = p.generate_folder_identifier(max_len=131)
    assert (
        len(folder_id) == 131
    ), f"Expected folder name length of 255 but got: {len(folder_id)}"


def test_generate_folder_identifier_max_len_exceeded():
    """Test folder identifier with max length exceeded."""
    p = Params({"key1": "value1" * 50, "key2": "value2" * 50})
    folder_id = p.generate_folder_identifier(max_len=50)
    assert folder_id.endswith(
        "_" + p.hash()[:8]
    ), f"Expected folder name to end with hash, but got: {folder_id[-8:]}"


def test_generate_folder_identifier_with_numeric_values():
    """Test folder identifier with numeric values in keys and values."""
    p = Params({"123": 456, "789": 101112})
    folder_id = p.generate_folder_identifier()
    expected_folder_id = "123-456_789-101112"
    assert (
        folder_id == expected_folder_id
    ), f"Expected '{expected_folder_id}' but got: {folder_id}"


def test_generate_folder_identifier_with_large_data():
    """Test folder identifier generation with a large number of params."""
    params = {f"key{i}": f"value{i}" for i in range(100)}
    p = Params(params)
    folder_id = p.generate_folder_identifier(max_len=1000)
    assert (
        len(folder_id) <= 1000
    ), f"Expected folder name length <= 1000 but got: {len(folder_id)}"


def test_generate_folder_identifier_empty_key_value():
    """Test folder identifier with empty key-value pairs."""
    p = Params({"key1": "", "key2": ""})
    folder_id = p.generate_folder_identifier()
    expected_folder_id = "key1-_key2-"
    assert (
        folder_id == expected_folder_id
    ), f"Expected '{expected_folder_id}' but got: {folder_id}"


def test_generate_folder_identifier_with_boolean_values():
    """Test folder identifier with boolean values."""
    p = Params({"key1": True, "key2": False})
    folder_id = p.generate_folder_identifier()
    expected_folder_id = "key1-True_key2-False"
    assert (
        folder_id == expected_folder_id
    ), f"Expected '{expected_folder_id}' but got: {folder_id}"
