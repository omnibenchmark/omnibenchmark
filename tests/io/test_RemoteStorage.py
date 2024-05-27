import pytest
from omni.io.RemoteStorage import RemoteStorage

def test_wrong_init():
    for ao in [1, 1.0, "1", [1], (1,), None]:
        with pytest.raises(ValueError):
            ss = RemoteStorage(ao, benchmark = "tb")
    
    ss = RemoteStorage({}, benchmark = "tb")
    assert type(ss) == RemoteStorage

    for ao in [1, 1.0, {}, [1], (1,), None, "", "-", ".", "0a"]:
        with pytest.raises(ValueError):
            ss = RemoteStorage({}, benchmark = ao)

    for ao in ["a","a0","A","A0","a_","a_0","A_","A_0"]:
        ss = RemoteStorage({}, benchmark = "a")
        assert type(ss) == RemoteStorage


def test_version():
    ss = RemoteStorage({}, benchmark = "a")

    ss.versions = ["0.1", "0.2", "1.0"]
    ss.set_current_version()
    assert ss.major_version == 1
    assert ss.minor_version == 0

    ss.set_current_version(major_version=1)
    assert ss.major_version == 1
    assert ss.minor_version == 0

    ss.set_current_version(major_version=0)
    assert ss.major_version == 0
    assert ss.minor_version == 2

    for mv in [2, "2", "1", ""]:
        with pytest.raises(ValueError):
            ss.set_current_version(major_version=mv)

    ss.set_current_version(minor_version=0)
    assert ss.major_version == 1
    assert ss.minor_version == 0

    for mv in [1, "1", "0", ""]:
        with pytest.raises(ValueError):
            ss.set_current_version(minor_version=mv)

    with pytest.raises(ValueError):
        ss.set_current_version(major_version=0, minor_version=0)

    ss.set_current_version(major_version=0, minor_version=1)
    assert ss.major_version == 0
    assert ss.minor_version == 1

    ss.set_current_version(major_version=0, minor_version=2)
    assert ss.major_version == 0
    assert ss.minor_version == 2

    ss.set_current_version(major_version=1, minor_version=0)
    assert ss.major_version == 1
    assert ss.minor_version == 0

def test_parse_version():
    ss = RemoteStorage({}, benchmark = "a")
    ss.versions = ["0.1", "0.2", "1.0"]
    ss.set_current_version()

    assert ss._parse_new_version(0, 1) == (0, 1)
    assert ss._parse_new_version(1, 1) == (1, 1)
    assert ss._parse_new_version(100, 100) == (100, 100)
    assert ss._parse_new_version(100, 0) == (100, 0)
    assert ss._parse_new_version(0, 100) == (0, 100)

    assert ss._parse_new_version(True, None) == (2, 0)
    assert ss._parse_new_version(True) == (2, 0)
    assert ss._parse_new_version(None, True) == (1, 1)
    assert ss._parse_new_version(True, True) == (2, 0)
    assert ss._parse_new_version(None, None) == (1, 1)

    for vv in [[1,"1"], ["1",1], [1,1.0], [1.0,1], [1.0,1.0], [1.0,"1.0"], ["1.0",1.0], ["1.0","1.0"],[None,0],[0,None],[False,None],[None,False],[False,False]]:
        with pytest.raises(ValueError):
            ss._parse_new_version(*vv)


def test_set_new_version():
    ss = RemoteStorage({}, benchmark = "a")
    ss.versions = ["0.1", "0.2", "1.0"]
    ss.set_current_version()

    ss.set_new_version()
    assert ss.major_version_new == 1 and ss.minor_version_new == 1

    ss.set_new_version(0,3)
    assert ss.major_version_new == 0 and ss.minor_version_new == 3

    ss.set_new_version(100,100)
    assert ss.major_version_new == 100 and ss.minor_version_new == 100

    ss.set_new_version(True)
    assert ss.major_version_new == 2 and ss.minor_version_new == 0

    ss.set_new_version(True, True)
    assert ss.major_version_new == 2 and ss.minor_version_new == 0

    ss.set_new_version(True, None)
    assert ss.major_version_new == 2 and ss.minor_version_new == 0

    ss.set_new_version(None, None)
    assert ss.major_version_new == 1 and ss.minor_version_new == 1

    ss.set_new_version(None, True)
    assert ss.major_version_new == 1 and ss.minor_version_new == 1

    # already exists
    with pytest.raises(ValueError):
        ss.set_new_version(1,0)




