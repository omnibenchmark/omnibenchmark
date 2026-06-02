"""Unit tests for filename truncation in benchmark._paths."""

from omnibenchmark.core._paths import (
    MAX_FILENAME_LEN,
    truncate_filename,
    truncate_path_filename,
)


def test_short_name_unchanged():
    assert truncate_filename("foo.txt") == "foo.txt"


def test_truncates_long_name_preserving_simple_extension():
    stem = "x" * 400
    name = stem + ".txt"
    out = truncate_filename(name)
    assert len(out) == MAX_FILENAME_LEN
    assert out.endswith(".txt")


def test_preserves_compound_tar_gz():
    stem = "dataset_" + "y" * 400
    name = stem + ".tar.gz"
    out = truncate_filename(name)
    assert len(out) == MAX_FILENAME_LEN
    assert out.endswith(".tar.gz")


def test_preserves_nii_gz():
    name = "scan_" + "z" * 300 + ".nii.gz"
    out = truncate_filename(name)
    assert len(out) == MAX_FILENAME_LEN
    assert out.endswith(".nii.gz")


def test_distinguishes_long_names_sharing_prefix():
    a = "shared_prefix_" + "a" * 300 + ".txt"
    b = "shared_prefix_" + "b" * 300 + ".txt"
    assert truncate_filename(a) != truncate_filename(b)


def test_dotted_stem_does_not_eat_filename():
    # All-suffix join would be ".v1.beta.release.candidate.txt" (>16 chars), so
    # the helper should fall back to just the trailing ".txt" extension.
    name = "my.v1.beta.release.candidate." + "x" * 300 + ".txt"
    out = truncate_filename(name)
    assert out.endswith(".txt")
    assert len(out) == MAX_FILENAME_LEN
    # The earlier dotted segments live in the stem, not the extension.
    assert "v1.beta" in out


def test_no_extension_still_truncates():
    name = "x" * 400
    out = truncate_filename(name)
    assert len(out) == MAX_FILENAME_LEN
    assert "." not in out  # no extension introduced


def test_truncate_path_filename_only_touches_basename():
    parent = "results/stage1/module/abcd1234"
    name = "x" * 400 + ".tar.gz"
    out = truncate_path_filename(f"{parent}/{name}")
    out_parent, _, out_name = out.rpartition("/")
    assert out_parent == parent
    assert out_name.endswith(".tar.gz")
    assert len(out_name) == MAX_FILENAME_LEN


def test_truncate_path_filename_short_path_unchanged():
    p = "results/stage1/module/abcd1234/foo.txt"
    assert truncate_path_filename(p) == p


def test_custom_limit():
    out = truncate_filename("a" * 50 + ".txt", limit=20)
    assert len(out) == 20
    assert out.endswith(".txt")
