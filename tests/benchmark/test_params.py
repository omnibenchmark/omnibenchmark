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
