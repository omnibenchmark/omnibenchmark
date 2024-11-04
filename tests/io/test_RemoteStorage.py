from unittest.mock import patch

import pytest
from packaging.version import InvalidVersion, Version

from omni.io.exception import RemoteStorageInvalidInputException
from omni.io.RemoteStorage import (
    DEFAULT_STORAGE_OPTIONS,
    RemoteStorage,
    is_valid_version,
)


class TestRemoteStorage:
    def test_is_valid_version(self):
        for ao in ["ab", [1], (1,), None, 1]:
            assert not is_valid_version(ao)
        for ao in ["1", "0.1"]:
            assert is_valid_version(ao)

    @patch.multiple(RemoteStorage, __abstractmethods__=set())
    def test__parse_benchmark(self):
        ss = RemoteStorage({}, benchmark="a")
        with pytest.raises(RemoteStorageInvalidInputException):
            for ao in [1, 0.1, [1], (1,), None]:
                ss._parse_benchmark(ao)
        assert ss.benchmark == "a"

    @patch.multiple(RemoteStorage, __abstractmethods__=set())
    def test__parse_auth_options(self):
        ss = RemoteStorage({}, benchmark="a")
        with pytest.raises(RemoteStorageInvalidInputException):
            for ao in [1, 0.1, [1], (1,), None]:
                ss._parse_auth_options(ao)
        assert ss.auth_options == {}

    @patch.multiple(RemoteStorage, __abstractmethods__=set())
    def test__parse_storage_options(self):
        ss = RemoteStorage({}, benchmark="a")
        with pytest.raises(RemoteStorageInvalidInputException):
            for ao in [1, 0.1, [1], (1,), None]:
                ss._parse_storage_options(ao)
        assert ss.storage_options == DEFAULT_STORAGE_OPTIONS

    @patch.multiple(RemoteStorage, __abstractmethods__=set())
    def test__parse_version(self):
        ss = RemoteStorage({}, benchmark="a")
        ss.versions = []
        assert ss._parse_version("0.1") == Version("0.1")
        assert ss._parse_version() == Version("0.1")
        assert ss._parse_version(None) == Version("0.1")
        assert ss._parse_version("0.2") == Version("0.2")

        ss.versions = [Version("0.1")]
        assert ss._parse_version() == Version("0.2")
        assert ss._parse_version(None) == Version("0.2")

    @patch.multiple(RemoteStorage, __abstractmethods__=set())
    def test_init_fails_with_invalid_arg_auth_options(self):
        for ao in [1, 1.0, "1", [1], (1,), None]:
            with pytest.raises(RemoteStorageInvalidInputException):
                ss = RemoteStorage(ao, benchmark="tb")

    @patch.multiple(RemoteStorage, __abstractmethods__=set())
    def test_init_success_with_empty_arg_auth_options(self):
        ss = RemoteStorage({}, benchmark="tb")
        assert type(ss) == RemoteStorage

    @patch.multiple(RemoteStorage, __abstractmethods__=set())
    def test_init_success_with_valid_arg_benchmark(self):
        for ao in ["a", "a0", "A", "A0", "a_", "a_0", "A_", "A_0"]:
            ss = RemoteStorage({}, benchmark="a")
            assert type(ss) == RemoteStorage

    @patch.multiple(RemoteStorage, __abstractmethods__=set())
    def test_set_version_success_with_valid_args(self):
        ss = RemoteStorage({}, benchmark="a")
        ss.versions = [Version("0.1"), Version("0.2"), Version("1.0")]
        ss.set_version()
        print(ss.version)
        assert ss.version == Version("1.1")

        ss = RemoteStorage({}, benchmark="a")
        ss.versions = [Version("0.1"), Version("0.2"), Version("1.0")]
        ss.set_version("1.0")
        print(ss.version)
        assert ss.version == Version("1.0")

        ss = RemoteStorage({}, benchmark="a")
        ss.versions = [Version("0.1"), Version("0.2"), Version("1.0")]
        print(ss.version)
        ss.set_version("0.1")
        assert ss.version == Version("0.1")
