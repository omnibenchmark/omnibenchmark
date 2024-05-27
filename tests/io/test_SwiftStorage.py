import pytest
from swiftclient.service import SwiftError, SwiftService
from omni.io.SwiftStorage import SwiftStorage 

class MockSwiftStorage():
    @staticmethod
    def connect():
        return MockSwiftServiceValid()

class MockSwiftServiceValid():
    def __enter__(self):
        return self
    def __exit__(self, exc_type, exc_value, traceback):
        pass

    def stat(self, container=None, objects=None):
        if container is None:
            return {'success': True}
        else:
            return [{'object': 'file1.txt', 'content_type': 'text/plain', 'hash': 'hash', 'last_modified': '2021-07-08T14:33:00.000000', 'bytes': 1, 'name': 'test_container', 'headers': {'x-object-meta-mtime': '1625731980.000000'}}, {'object': 'file2.txt', 'content_type': 'text/plain', 'hash': 'hash', 'last_modified': '2021-07-08T14:33:00.000000', 'bytes': 1, 'name': 'test_container', 'headers': {'x-object-meta-mtime': '1625731980.000000'}}]

    def list(self, container=None, options=None):
        if container is None:
            return [{'listing': [{'name': 'test.0.1', 'count': 1, 'bytes': 1, 'last_modified': '2021-07-08T14:33:00.000000', 'hash': 'hash'},{'name': 'test.0.2', 'count': 1, 'bytes': 1, 'last_modified': '2021-07-08T14:33:00.000000', 'hash': 'hash'},{'name': 'test.1.0', 'count': 1, 'bytes': 1, 'last_modified': '2021-07-08T14:33:00.000000', 'hash': 'hash'},{'name': 'test2.0.1', 'count': 1, 'bytes': 1, 'last_modified': '2021-07-08T14:33:00.000000', 'hash': 'hash'}]}]
        else:
            return [{'listing': [{'name': 'file1.txt', 'count': 1, 'bytes': 1, 'last_modified': '2021-07-08T14:33:00.000000', 'hash': 'hash'}, {'name': 'file2.txt', 'count': 1, 'bytes': 1, 'last_modified': '2021-07-08T14:33:00.000000', 'hash': 'hash'}]}]
    
    def post(self, container, options = None):
        return {'success': True, 'container': 'tb.0.1', 'object': None, 'headers': {'X-Container-Read': '.r:*,.rlistings'}, 'action': 'post_container'}

    def copy(self,container, objects, destination):
        return [{'success': True, 'action': 'copy_object', 'object': 'file1.txt', 'destination': 'test.1.1'}, {'success': True, 'action': 'copy_object', 'object': 'file2.txt', 'destination': 'test.1.1'}]

    def delete(self,container, objects):
        return [{'success': True, 'action': 'copy_object', 'object': 'file1.txt', 'destination': 'test.1.1'}, {'success': True, 'action': 'copy_object', 'object': 'file2.txt', 'destination': 'test.1.1'}]


def mock_test_connect(*args, **kwargs):
    return MockSwiftServiceValid()

def test__test_connect(monkeypatch):
    monkeypatch.setattr(SwiftStorage, "connect", mock_test_connect)

    ss = SwiftStorage(auth_options={'preauthurl': "https://", "preauthtoken": "token"}, benchmark="test")
    assert ss.benchmark == "test"
    ss._test_connect()

def test_init(monkeypatch):
    monkeypatch.setattr(SwiftStorage, "connect", mock_test_connect)

    ss = SwiftStorage(auth_options={'preauthurl': "https://", "preauthtoken": "token"}, benchmark="test")

def test_init_public():
    ss = SwiftStorage(auth_options={'preauthurl': 'https://cloud.s3it.uzh.ch:8080/v1/AUTH_fed6cdab6446481597fad8aa67751bd5', 'retries': 3, 'starting_backoff': 1, 'max_backoff': 8}, benchmark="tb2")


def test__get_containers(monkeypatch):
    monkeypatch.setattr(SwiftStorage, "connect", mock_test_connect)

    ss = SwiftStorage(auth_options={'preauthurl': "https://", "preauthtoken": "token"}, benchmark="test")
    ss._get_containers()
    assert ss.containers == ['test.0.1', 'test.0.2', 'test.1.0', 'test2.0.1']
    
def test__get_benchmarks(monkeypatch):
    monkeypatch.setattr(SwiftStorage, "connect", mock_test_connect)

    ss = SwiftStorage(auth_options={'preauthurl': "https://", "preauthtoken": "token"}, benchmark="test")
    ss._get_benchmarks()
    assert ss.benchmarks == ['test', 'test2']

    ss._get_benchmarks(update=False)
    assert ss.benchmarks == ['test', 'test2']


def test__create_benchmark(monkeypatch):
    monkeypatch.setattr(SwiftStorage, "connect", mock_test_connect)

    ss = SwiftStorage(auth_options={'preauthurl': "https://", "preauthtoken": "token"}, benchmark="test")
    ss._create_benchmark('test3')
    # assert ss.benchmarks == ['test', 'test2','test3']

    ss = SwiftStorage(auth_options={'preauthurl': "https://", "preauthtoken": "token"}, benchmark="test")
    ss._create_benchmark('test3', update=False)


def test__get_versions(monkeypatch):
    monkeypatch.setattr(SwiftStorage, "connect", mock_test_connect)

    ss = SwiftStorage(auth_options={'preauthurl': "https://", "preauthtoken": "token"}, benchmark="test")
    ss._get_versions()
    assert ss.versions == ['0.1', '0.2', '1.0']

    ss._get_versions(update=False)
    assert ss.versions == ['0.1', '0.2', '1.0']

def test__get_versions_public():
    ss = SwiftStorage(auth_options={'preauthurl': 'https://cloud.s3it.uzh.ch:8080/v1/AUTH_fed6cdab6446481597fad8aa67751bd5', 'retries': 3, 'starting_backoff': 1, 'max_backoff': 8}, benchmark="tb2")
    ss._get_versions()
    assert ss.versions == ['0.1']

    ss._get_versions(readonly=True)
    assert ss.versions == ['0.1']

def test__create_new_version(monkeypatch):
    monkeypatch.setattr(SwiftStorage, "connect", mock_test_connect)

    ss = SwiftStorage(auth_options={'preauthurl': "https://", "preauthtoken": "token"}, benchmark="test")
    ss.set_new_version()
    assert ss.major_version == 1
    assert ss.minor_version == 0
    assert ss.major_version_new == 1
    assert ss.minor_version_new == 1
    assert ss.versions == ['0.1', '0.2', '1.0']

    # fails because no new version is created in mock test
    with pytest.raises(ValueError):
        ss._create_new_version()


def test__set_current_version(monkeypatch):
    monkeypatch.setattr(SwiftStorage, "connect", mock_test_connect)

    ss = SwiftStorage(auth_options={'preauthurl': "https://", "preauthtoken": "token"}, benchmark="test")
    ss.set_current_version()
    assert ss.major_version == 1
    assert ss.minor_version == 0

def test__set_current_version_public():
    ss = SwiftStorage(auth_options={'preauthurl': 'https://cloud.s3it.uzh.ch:8080/v1/AUTH_fed6cdab6446481597fad8aa67751bd5', 'retries': 3, 'starting_backoff': 1, 'max_backoff': 8}, benchmark="tb2")
    ss.set_current_version()
    assert ss.major_version == 0
    assert ss.minor_version == 1

def test__set_new_version(monkeypatch):
    monkeypatch.setattr(SwiftStorage, "connect", mock_test_connect)

    ss = SwiftStorage(auth_options={'preauthurl': "https://", "preauthtoken": "token"}, benchmark="test")
    ss.set_new_version()
    assert ss.major_version_new == 1
    assert ss.minor_version_new == 1


def test__get_objects(monkeypatch):
    monkeypatch.setattr(SwiftStorage, "connect", mock_test_connect)

    ss = SwiftStorage(auth_options={'preauthurl': "https://", "preauthtoken": "token"}, benchmark="test")
    ss.set_current_version()
    ss._get_objects()
    assert ss.files.keys() == {'file1.txt', 'file2.txt'}
    assert ss.files['file1.txt'].keys() == {'hash', 'last_modified', 'size', 'x-object-meta-mtime', 'symlink_path'}
    assert ss.files['file2.txt'].keys() == {'hash', 'last_modified', 'size', 'x-object-meta-mtime', 'symlink_path'}
    assert ss.files['file1.txt']['size'] == 1
    assert ss.files['file2.txt']['size'] == 1
    assert ss.files['file1.txt']['x-object-meta-mtime'] == '1625731980.000000'
    assert ss.files['file2.txt']['x-object-meta-mtime'] == '1625731980.000000'
    assert ss.files['file1.txt']['hash'] == 'hash'
    assert ss.files['file2.txt']['hash'] == 'hash'
    assert ss.files['file1.txt']['last_modified'] == '2021-07-08T14:33:00.000000'
    assert ss.files['file2.txt']['last_modified'] == '2021-07-08T14:33:00.000000'

def test__get_objects_public():
    ss = SwiftStorage(auth_options={'preauthurl': 'https://cloud.s3it.uzh.ch:8080/v1/AUTH_fed6cdab6446481597fad8aa67751bd5', 'retries': 3, 'starting_backoff': 1, 'max_backoff': 8}, benchmark="tb2")
    with pytest.raises(ValueError):
        ss._get_objects()

    ss.set_current_version()
    ss._get_objects()
    assert ss.files.keys() == {'test.txt', 'test2.txt'}
    assert ss.files['test.txt'].keys() == {'hash', 'last_modified', 'size', 'x-object-meta-mtime', 'symlink_path', 'accesstime'}


def test_find_objects_to_copy(monkeypatch):
    monkeypatch.setattr(SwiftStorage, "connect", mock_test_connect)

    ss = SwiftStorage(auth_options={'preauthurl': "https://", "preauthtoken": "token"}, benchmark="test")
    ss.set_current_version()
    ss.find_objects_to_copy(tagging_type="all")
    assert ss.files.keys() == {'file1.txt', 'file2.txt'}
    assert ss.files['file1.txt'].keys() == {'hash', 'last_modified', 'size', 'x-object-meta-mtime', 'symlink_path', 'copy'}
    assert ss.files['file2.txt'].keys() == {'hash', 'last_modified', 'size', 'x-object-meta-mtime', 'symlink_path', 'copy'}
    assert type(ss.files['file1.txt']['copy']) == bool
    assert type(ss.files['file2.txt']['copy']) == bool

    with pytest.raises(ValueError):
        ss.find_objects_to_copy(tagging_type="other")

    with pytest.raises(ValueError):
        ss.find_objects_to_copy(reference_time="2024")
    
    import datetime
    ss.find_objects_to_copy(reference_time=datetime.datetime.now())
    assert ss.files.keys() == {'file1.txt', 'file2.txt'}
    assert ss.files['file1.txt'].keys() == {'hash', 'last_modified', 'size', 'x-object-meta-mtime', 'symlink_path', 'copy'}
    assert ss.files['file2.txt'].keys() == {'hash', 'last_modified', 'size', 'x-object-meta-mtime', 'symlink_path', 'copy'}
    assert ss.files['file1.txt']['copy'] == True
    assert ss.files['file2.txt']['copy'] == True



def test_copy_objects(monkeypatch):
    monkeypatch.setattr(SwiftStorage, "connect", mock_test_connect)

    ss = SwiftStorage(auth_options={'preauthurl': "https://", "preauthtoken": "token"}, benchmark="test")
    ss.set_new_version()
    ss.find_objects_to_copy()
    ss.copy_objects()
    assert ss.files.keys() == {'file1.txt', 'file2.txt'}
    assert ss.files['file1.txt'].keys() == {'hash', 'last_modified', 'size', 'x-object-meta-mtime', 'symlink_path', 'copy', 'copied'}
    assert ss.files['file2.txt'].keys() == {'hash', 'last_modified', 'size', 'x-object-meta-mtime', 'symlink_path', 'copy', 'copied'}
    assert type(ss.files['file1.txt']['copied']) == bool
    assert type(ss.files['file2.txt']['copied']) == bool

    with pytest.raises(NotImplementedError):
        ss.copy_objects(type="symlink")
    with pytest.raises(ValueError):
        ss.copy_objects(type="other")



def test_create_new_version(monkeypatch):
    monkeypatch.setattr(SwiftStorage, "connect", mock_test_connect)

    ss = SwiftStorage(auth_options={'preauthurl': "https://", "preauthtoken": "token"}, benchmark="test")

    # fails because no new version is created in mock test
    with pytest.raises(ValueError):
        ss.create_new_version()
