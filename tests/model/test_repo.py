import pytest
import hashlib

from omnibenchmark.model.repo import get_repo_hash


@pytest.mark.short
def test_get_repo_hash_basic():
    """Test basic functionality with valid inputs."""
    repo_url = "https://github.com/user/repo.git"
    commit_hash = "abc123def456"

    result = get_repo_hash(repo_url, commit_hash)

    # Should return a SHA256 hash (64 characters, hexadecimal)
    assert len(result) == 64
    assert all(c in "0123456789abcdef" for c in result)

    # Verify it's the correct hash
    expected = "d8af15088c0b9c961efcadeda3faaf0a07cda14c2bba079b4e69bbb46e2f4066"
    assert result == expected


@pytest.mark.short
def test_get_repo_hash_consistency():
    """Test that same inputs always produce same hash."""
    repo_url = "https://github.com/user/repo.git"
    commit_hash = "abc123def456"

    hash1 = get_repo_hash(repo_url, commit_hash)
    hash2 = get_repo_hash(repo_url, commit_hash)

    assert hash1 == hash2


@pytest.mark.short
def test_get_repo_hash_different_inputs():
    """Test that different inputs produce different hashes."""
    repo_url1 = "https://github.com/user/repo1.git"
    repo_url2 = "https://github.com/user/repo2.git"
    commit_hash = "abc123def456"

    hash1 = get_repo_hash(repo_url1, commit_hash)
    hash2 = get_repo_hash(repo_url2, commit_hash)

    assert hash1 != hash2


@pytest.mark.short
def test_get_repo_hash_different_commits():
    """Test that different commits produce different hashes."""
    repo_url = "https://github.com/user/repo.git"
    commit_hash1 = "abc123def456"
    commit_hash2 = "def456abc123"

    hash1 = get_repo_hash(repo_url, commit_hash1)
    hash2 = get_repo_hash(repo_url, commit_hash2)

    assert hash1 != hash2


@pytest.mark.short
def test_get_repo_hash_none_repo_url():
    """Test that None repo_url raises ValueError."""
    with pytest.raises(
        ValueError, match="Both repo_url and commit_hash must be provided"
    ):
        get_repo_hash(None, "abc123")


@pytest.mark.short
def test_get_repo_hash_none_commit_hash():
    """Test that None commit_hash raises ValueError."""
    with pytest.raises(
        ValueError, match="Both repo_url and commit_hash must be provided"
    ):
        get_repo_hash("https://github.com/user/repo.git", None)


@pytest.mark.short
def test_get_repo_hash_both_none():
    """Test that both None values raise ValueError."""
    with pytest.raises(
        ValueError, match="Both repo_url and commit_hash must be provided"
    ):
        get_repo_hash(None, None)


@pytest.mark.short
def test_get_repo_hash_empty_repo_url():
    """Test that empty repo_url raises ValueError."""
    with pytest.raises(
        ValueError, match="Both repo_url and commit_hash must be provided"
    ):
        get_repo_hash("", "abc123")


@pytest.mark.short
def test_get_repo_hash_empty_commit_hash():
    """Test that empty commit_hash raises ValueError."""
    with pytest.raises(
        ValueError, match="Both repo_url and commit_hash must be provided"
    ):
        get_repo_hash("https://github.com/user/repo.git", "")


@pytest.mark.short
def test_get_repo_hash_both_empty():
    """Test that both empty values raise ValueError."""
    with pytest.raises(
        ValueError, match="Both repo_url and commit_hash must be provided"
    ):
        get_repo_hash("", "")


@pytest.mark.short
def test_get_repo_hash_various_url_formats():
    """Test with different repository URL formats."""
    commit_hash = "abc123def456"

    urls = [
        "https://github.com/user/repo.git",
        "https://github.com/user/repo",
        "git@github.com:user/repo.git",
        "https://gitlab.com/user/repo.git",
        "https://bitbucket.org/user/repo.git",
        "file:///path/to/local/repo",
    ]

    hashes = []
    for url in urls:
        hash_value = get_repo_hash(url, commit_hash)
        assert len(hash_value) == 64
        assert all(c in "0123456789abcdef" for c in hash_value)
        hashes.append(hash_value)

    # All hashes should be different
    assert len(set(hashes)) == len(hashes)


@pytest.mark.short
def test_get_repo_hash_various_commit_formats():
    """Test with different commit hash formats."""
    repo_url = "https://github.com/user/repo.git"

    commits = [
        "abc123def456",  # Short hash
        "abc123def456789012345678901234567890abcd",  # Full SHA-1 hash
        "v1.0.0.final",  # Tag (8+ chars)
        "main-branch",  # Branch name (8+ chars)
        "feature/new-feature",  # Branch with slash
        "release-1.2.3",  # Branch with dashes
    ]

    hashes = []
    for commit in commits:
        hash_value = get_repo_hash(repo_url, commit)
        assert len(hash_value) == 64
        assert all(c in "0123456789abcdef" for c in hash_value)
        hashes.append(hash_value)

    # All hashes should be different
    assert len(set(hashes)) == len(hashes)


@pytest.mark.short
def test_get_repo_hash_special_characters():
    """Test handling of special characters in URLs and commits."""
    # URL with query parameters and special chars
    repo_url = "https://github.com/user/repo.git?ref=main&token=abc123"
    commit_hash = "feature/update-deps#123"

    result = get_repo_hash(repo_url, commit_hash)

    assert len(result) == 64
    assert all(c in "0123456789abcdef" for c in result)

    # Verify expected hash calculation
    expected = hashlib.sha256(f"{repo_url}@{commit_hash}".encode()).hexdigest()
    assert result == expected


@pytest.mark.short
def test_get_repo_hash_unicode_characters():
    """Test handling of unicode characters."""
    repo_url = "https://github.com/用户/repository.git"
    commit_hash = "commit-with-émoji-🚀"

    result = get_repo_hash(repo_url, commit_hash)

    assert len(result) == 64
    assert all(c in "0123456789abcdef" for c in result)

    # Verify expected hash calculation
    expected = hashlib.sha256(f"{repo_url}@{commit_hash}".encode()).hexdigest()
    assert result == expected


@pytest.mark.short
def test_get_repo_hash_whitespace_handling():
    """Test that whitespace is preserved in hash calculation."""
    repo_url1 = "https://github.com/user/repo.git"
    repo_url2 = " https://github.com/user/repo.git "
    commit_hash = "abc123def456"

    hash1 = get_repo_hash(repo_url1, commit_hash)
    hash2 = get_repo_hash(repo_url2, commit_hash)

    # Whitespace should make hashes different
    assert hash1 != hash2


@pytest.mark.short
def test_get_repo_hash_very_long_inputs():
    """Test with very long repository URLs and commit hashes."""
    # Very long repository URL
    long_repo_url = "https://github.com/" + "a" * 1000 + "/repo.git"
    # Very long commit hash
    long_commit_hash = "b" * 1000

    result = get_repo_hash(long_repo_url, long_commit_hash)

    assert len(result) == 64
    assert all(c in "0123456789abcdef" for c in result)

    # Verify expected hash calculation
    expected = hashlib.sha256(
        f"{long_repo_url}@{long_commit_hash}".encode()
    ).hexdigest()
    assert result == expected


@pytest.mark.short
def test_get_repo_hash_format_string():
    """Test the specific format used for hash generation."""
    repo_url = "https://github.com/user/repo.git"
    commit_hash = "abc123def456"

    # The function should hash "repo_url@commit_hash"
    expected_input = f"{repo_url}@{commit_hash}"
    expected_hash = hashlib.sha256(expected_input.encode()).hexdigest()

    result = get_repo_hash(repo_url, commit_hash)

    assert result == expected_hash


@pytest.mark.short
def test_get_repo_hash_short_commit_hash():
    """Test that short commit hashes work correctly (matching original utils.py behavior)."""
    repo_url = "https://github.com/user/repo.git"

    short_commits = [
        "a",
        "ab",
        "abc",
        "abcd",
        "abcde",
        "abcdef",
        "abcdefg",
    ]

    for commit in short_commits:
        result = get_repo_hash(repo_url, commit)
        assert len(result) == 64
        assert all(c in "0123456789abcdef" for c in result)
        # Verify expected hash calculation
        expected = hashlib.sha256(f"{repo_url}@{commit}".encode()).hexdigest()
        assert result == expected


@pytest.mark.short
def test_get_repo_hash_various_commit_lengths():
    """Test that commits of various lengths work correctly."""
    repo_url = "https://github.com/user/repo.git"

    commits = [
        "a",  # 1 char
        "dev",  # 3 chars
        "v1.0",  # 4 chars
        "main",  # 4 chars
        "abcdefg",  # 7 chars
        "abcdefgh",  # 8 chars
        "123456789",  # 9 chars
    ]

    hashes = []
    for commit in commits:
        result = get_repo_hash(repo_url, commit)
        assert len(result) == 64
        assert all(c in "0123456789abcdef" for c in result)
        hashes.append(result)

    # All hashes should be different
    assert len(set(hashes)) == len(hashes)
