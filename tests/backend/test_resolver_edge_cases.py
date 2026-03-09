"""
Edge case and corner case unit tests for ModuleResolver.

Tests focus on boundary conditions and unusual input combinations
through the public API only. Tests are marked as 'short' since they
test logic without heavy I/O.
"""

import pytest
from pathlib import Path

from omnibenchmark.backend.resolver import ModuleResolver
from omnibenchmark.model import (
    Module,
    Repository,
    SoftwareBackendEnum,
)


# ---------------------------------------------------------------------------
# Test Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def minimal_module():
    """Minimal Module for testing."""
    return Module(
        id="test_module",
        repository=Repository(
            url="https://github.com/test/repo.git",
            commit="abc123",
        ),
        software_environment="test_env",
    )


@pytest.fixture
def resolver(tmp_path):
    """Basic ModuleResolver instance."""
    return ModuleResolver(
        work_base_dir=tmp_path / "work",
        cache_dir=tmp_path / "cache",
        output_dir=tmp_path / "output",
    )


# ---------------------------------------------------------------------------
# Test ModuleResolver initialization edge cases
# ---------------------------------------------------------------------------


@pytest.mark.short
class TestModuleResolverInitialization:
    """Test ModuleResolver initialization edge cases."""

    def test_default_directories_used(self):
        """Test that default directories are used when not specified."""
        resolver = ModuleResolver()

        assert resolver.work_base_dir == Path(".modules")
        assert resolver.cache_dir is None
        assert resolver.output_dir == Path(".")

    def test_custom_directories(self, tmp_path):
        """Test custom directory configuration."""
        work_dir = tmp_path / "custom_work"
        cache_dir = tmp_path / "custom_cache"
        output_dir = tmp_path / "custom_output"

        resolver = ModuleResolver(
            work_base_dir=work_dir,
            cache_dir=cache_dir,
            output_dir=output_dir,
        )

        assert resolver.work_base_dir == work_dir
        assert resolver.cache_dir == cache_dir
        assert resolver.output_dir == output_dir

    def test_software_environments_empty_dict_default(self):
        """Test that software_environments defaults to empty dict."""
        resolver = ModuleResolver()

        assert resolver.software_environments == {}

    def test_software_environments_empty_dict_allowed(self):
        """Test providing empty software environments dict."""
        envs = {}

        resolver = ModuleResolver(software_environments=envs)

        assert resolver.software_environments == envs

    def test_dirty_copy_cache_initialized_empty(self):
        """Test that dirty copy cache is initialized empty."""
        resolver = ModuleResolver()

        assert resolver._dirty_copy_cache == {}

    def test_benchmark_dir_default(self):
        """Test that benchmark_dir defaults to current directory."""
        resolver = ModuleResolver()

        assert resolver.benchmark_dir == Path(".")

    def test_software_backend_none_default(self):
        """Test that software_backend defaults to None."""
        resolver = ModuleResolver()

        assert resolver.software_backend is None

    def test_multiple_software_backends_supported(self):
        """Test that all software backend types can be configured."""
        for backend in SoftwareBackendEnum:
            resolver = ModuleResolver(software_backend=backend)
            assert resolver.software_backend == backend


# ---------------------------------------------------------------------------
# Test path handling edge cases
# ---------------------------------------------------------------------------


@pytest.mark.short
class TestPathHandlingEdgeCases:
    """Test edge cases in path handling."""

    def test_work_dir_with_spaces(self, tmp_path):
        """Test work directory path with spaces."""
        work_dir = tmp_path / "work dir with spaces"
        resolver = ModuleResolver(work_base_dir=work_dir)

        assert resolver.work_base_dir == work_dir

    def test_work_dir_with_unicode(self, tmp_path):
        """Test work directory path with Unicode characters."""
        work_dir = tmp_path / "work_测试_ñ"
        resolver = ModuleResolver(work_base_dir=work_dir)

        assert resolver.work_base_dir == work_dir

    def test_very_deep_directory_structure(self, tmp_path):
        """Test very deep directory nesting."""
        deep_path = tmp_path
        for i in range(20):
            deep_path = deep_path / f"level{i}"

        resolver = ModuleResolver(work_base_dir=deep_path)
        assert resolver.work_base_dir == deep_path

    def test_relative_work_dir(self):
        """Test relative work directory path."""
        resolver = ModuleResolver(work_base_dir=Path("./relative/path"))

        assert resolver.work_base_dir == Path("./relative/path")

    def test_absolute_work_dir(self, tmp_path):
        """Test absolute work directory path."""
        abs_path = tmp_path.absolute()
        resolver = ModuleResolver(work_base_dir=abs_path)

        assert resolver.work_base_dir == abs_path


# ---------------------------------------------------------------------------
# Test module configuration edge cases
# ---------------------------------------------------------------------------


@pytest.mark.short
class TestModuleConfigurationEdgeCases:
    """Test edge cases in module configuration."""

    def test_module_with_very_long_id(self):
        """Test module with extremely long ID."""
        long_id = "a" * 500
        module = Module(
            id=long_id,
            repository=Repository(
                url="https://github.com/test/repo.git",
                commit="abc123",
            ),
            software_environment="test",
        )

        # Should create without error
        assert module.id == long_id

    def test_module_with_special_chars_in_id(self):
        """Test module with special characters in ID."""
        module = Module(
            id="test-module@v1.0",
            repository=Repository(
                url="https://github.com/test/repo.git",
                commit="abc123",
            ),
            software_environment="test",
        )

        # Should create without error
        assert module.id == "test-module@v1.0"

    def test_module_with_unicode_in_id(self):
        """Test module with Unicode in ID."""
        module = Module(
            id="test_模块_αβγ",
            repository=Repository(
                url="https://github.com/test/repo.git",
                commit="abc123",
            ),
            software_environment="test",
        )

        # Should create without error
        assert "模块" in module.id

    def test_module_with_very_long_commit_hash(self):
        """Test module with very long commit hash."""
        long_hash = "a" * 100
        module = Module(
            id="test",
            repository=Repository(
                url="https://github.com/test/repo.git",
                commit=long_hash,
            ),
            software_environment="test",
        )

        assert module.repository.commit == long_hash

    def test_module_with_various_url_formats(self):
        """Test module with different URL formats."""
        url_formats = [
            "https://github.com/org/repo.git",
            "git@github.com:org/repo.git",
            "ssh://git@github.com/org/repo.git",
            "https://gitlab.com/group/subgroup/repo.git",
            "file:///local/path/to/repo",
        ]

        for url in url_formats:
            module = Module(
                id="test",
                repository=Repository(url=url, commit="abc123"),
                software_environment="test",
            )
            assert module.repository.url == url


# ---------------------------------------------------------------------------
# Test software backend configuration edge cases
# ---------------------------------------------------------------------------


@pytest.mark.short
class TestSoftwareBackendConfigurationEdgeCases:
    """Test edge cases in software backend configuration."""

    def test_resolver_with_all_backend_types(self, tmp_path):
        """Test that resolver accepts all backend types."""
        for backend in SoftwareBackendEnum:
            resolver = ModuleResolver(
                work_base_dir=tmp_path,
                software_backend=backend,
            )
            assert resolver.software_backend == backend

    def test_resolver_without_software_backend(self, tmp_path):
        """Test resolver without specifying software backend."""
        resolver = ModuleResolver(work_base_dir=tmp_path)
        assert resolver.software_backend is None


# ---------------------------------------------------------------------------
# Test dirty copy cache edge cases
# ---------------------------------------------------------------------------


@pytest.mark.short
class TestDirtyCopyCacheEdgeCases:
    """Test edge cases in dirty copy caching."""

    def test_cache_is_per_resolver_instance(self, tmp_path):
        """Test that dirty copy cache is per-resolver instance."""
        resolver1 = ModuleResolver(work_base_dir=tmp_path / "r1")
        resolver2 = ModuleResolver(work_base_dir=tmp_path / "r2")

        # Caches should be independent
        assert resolver1._dirty_copy_cache is not resolver2._dirty_copy_cache

    def test_cache_starts_empty(self, tmp_path):
        """Test that cache starts empty for each resolver."""
        resolver = ModuleResolver(work_base_dir=tmp_path)

        # Cache should be empty initially
        assert len(resolver._dirty_copy_cache) == 0

    def test_multiple_resolvers_independent_caches(self, tmp_path):
        """Test that multiple resolvers maintain independent caches."""
        resolvers = [ModuleResolver(work_base_dir=tmp_path / f"r{i}") for i in range(5)]

        # All caches should be independent
        cache_ids = [id(r._dirty_copy_cache) for r in resolvers]
        assert len(set(cache_ids)) == 5
