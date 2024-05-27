import pytest
import re
import requests
from omni.io.S3Storage import S3Storage


class MockS3Storage:
    def connect():
        return MockS3ServiceValid()


class MockbucketsCollectionManager:
    def __init__(self):
        self.buckets = [
            MockS3Bucket("test.0.1"),
            MockS3Bucket("test.0.2"),
            MockS3Bucket("test.1.0"),
            MockS3Bucket("test2.0.1"),
        ]

    def all(self):
        return self.buckets

    def iterator(self):
        return self.buckets


class MockS3Bucket:
    def __init__(self, name):
        self.name = name
        self.creation_date = "2021-07-08T14:33:00.000000"
        self.objects = MockS3objectsCollectionManager()

    def create_bucket(self):
        return {"success": True, "bucket": self.name, "action": "create_bucket"}

    def put_object(self, Key, Body):
        return {
            "success": True,
            "bucket": self.name,
            "object": Key,
            "action": "put_object",
        }

    def delete_objects(self, Delete):
        return {
            "success": True,
            "bucket": self.name,
            "objects": Delete,
            "action": "delete_objects",
        }

    def Policy(self):
        return MockS3Policy(self.name)


class MockS3Policy:
    def __init__(self, bucket):
        self.bucket = bucket

    def put(self, Policy):
        return {"success": True, "bucket": self.bucket, "action": "put_bucket_policy"}

    def delete(self):
        return {
            "success": True,
            "bucket": self.bucket,
            "action": "delete_bucket_policy",
        }


class MockS3objectsCollectionManager:
    def __init__(self):
        self.objects = [MockS3Object("file1.txt"), MockS3Object("file2.txt")]

    def all(self):
        return self.objects

    def iterator(self):
        return self.objects


class MockS3Object:
    def __init__(self, key):
        self.key = key
        self.last_modified = "2021-07-08T14:33:00.000000"
        self.size = 1
        self.e_tag = "hash"
        self.metadata = {"x-object-meta-mtime": "1625731980.000000"}

    def copy_from(self, CopySource):
        return {"success": True, "object": self.key, "action": "copy_object"}


class MockS3ServiceValid:
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        pass

    def __init__(self):
        self.buckets = MockbucketsCollectionManager()

    def create_bucket(self, Bucket):
        return {"success": True, "bucket": Bucket, "action": "create_bucket"}

    def Bucket(self, name):
        return MockS3Bucket(name)

    def Object(self, bucket, key):
        return MockS3Object(key)


def mock_test_connect(*args, **kwargs):
    return MockS3ServiceValid()


class MockResponse:
    def __init__(self, url):
        self.status_code = 200
        self.ok = True

        self.headers = {
            "Accept-Ranges": "bytes",
            "Content-Length": "2",
            "Content-Type": "application/octet-stream",
            "ETag": '"49f68a5c8493ec2c0bf489821c21fc3b"',
            "Last-Modified": "Fri, 03 May 2024 10:44:26 GMT",
            "Server": "MinIO",
            "Strict-Transport-Security": "max-age=31536000; includeSubDomains",
            "Vary": "Origin, Accept-Encoding",
            "X-Amz-Id-2": "dd9025bab4ad464b049177c95eb6ebf374d3b3fd1af9251148b658df7ac2e3e8",
            "X-Amz-Request-Id": "17D045E2CA527BB1",
            "X-Content-Type-Options": "nosniff",
            "X-Xss-Protection": "1; mode=block",
            "x-amz-tagging-count": "2",
            "Date": "Fri, 17 May 2024 12:18:18 GMT",
        }

        if re.search(f"overview", url):
            self.text = '<?xml version="1.0" encoding="UTF-8"?>\n<ListBucketResult xmlns="http://s3.amazonaws.com/doc/2006-03-01/"><Name>tb.overview</Name><Prefix></Prefix><Marker></Marker><MaxKeys>1000</MaxKeys><IsTruncated>false</IsTruncated><Contents><Key>0.1</Key><LastModified>2024-05-06T10:11:14.312Z</LastModified><ETag>&#34;d41d8cd98f00b204e9800998ecf8427e&#34;</ETag><Size>0</Size><Owner><ID>02d6176db174dc93cb1b899f7c6078f08654445fe8cf1b6ce98d8855f66bdbf4</ID><DisplayName>minio</DisplayName></Owner><StorageClass>STANDARD</StorageClass></Contents><Contents><Key>0.2</Key><LastModified>2024-05-06T11:29:30.299Z</LastModified><ETag>&#34;d41d8cd98f00b204e9800998ecf8427e&#34;</ETag><Size>0</Size><Owner><ID>02d6176db174dc93cb1b899f7c6078f08654445fe8cf1b6ce98d8855f66bdbf4</ID><DisplayName>minio</DisplayName></Owner><StorageClass>STANDARD</StorageClass></Contents></ListBucketResult>'
        else:  # re.match(f"overview", url):
            self.text = '<?xml version="1.0" encoding="UTF-8"?>\n<ListBucketResult xmlns="http://s3.amazonaws.com/doc/2006-03-01/"><Name>tb.0.1</Name><Prefix></Prefix><Marker></Marker><MaxKeys>1000</MaxKeys><IsTruncated>false</IsTruncated><Contents><Key>test.txt</Key><LastModified>2024-05-02T14:46:22.966Z</LastModified><ETag>&#34;d636b8ffe4e3b6ee520fa65691124a10&#34;</ETag><Size>338</Size><Owner><ID>02d6176db174dc93cb1b899f7c6078f08654445fe8cf1b6ce98d8855f66bdbf4</ID><DisplayName>minio</DisplayName></Owner><StorageClass>STANDARD</StorageClass></Contents><Contents><Key>test2.txt</Key><LastModified>2024-05-02T14:46:22.893Z</LastModified><ETag>&#34;3aafbd94758e7b9fb6e7cb1d295c1c22&#34;</ETag><Size>338</Size><Owner><ID>02d6176db174dc93cb1b899f7c6078f08654445fe8cf1b6ce98d8855f66bdbf4</ID><DisplayName>minio</DisplayName></Owner><StorageClass>STANDARD</StorageClass></Contents></ListBucketResult>'


def mock_test_requests_get(url, **kwargs):
    return MockResponse(url)


def test__test_connect(monkeypatch):
    monkeypatch.setattr(S3Storage, "connect", mock_test_connect)

    ss = S3Storage(
        auth_options={
            "endpoint": "http://omnibenchmark.mls.uzh.ch:9000",
            "access_key": "asdfasdf",
            "secret_key": "asdfasdfasdf",
            "secure": False,
        },
        benchmark="test",
    )
    assert ss.benchmark == "test"
    ss._test_connect()


def test_init(monkeypatch):
    monkeypatch.setattr(S3Storage, "connect", mock_test_connect)

    ss = S3Storage(
        auth_options={
            "endpoint": "http://omnibenchmark.mls.uzh.ch:9000",
            "access_key": "asdfasdf",
            "secret_key": "asdfasdfasdf",
            "secure": False,
        },
        benchmark="test",
    )


def test_init_public(monkeypatch):
    monkeypatch.setattr(requests, "get", mock_test_requests_get)

    ss = S3Storage(
        auth_options={
            "endpoint": "http://omnibenchmark.mls.uzh.ch:9000",
            "secure": False,
        },
        benchmark="test",
    )


def test__get_containers(monkeypatch):
    monkeypatch.setattr(S3Storage, "connect", mock_test_connect)

    ss = S3Storage(
        auth_options={
            "endpoint": "http://omnibenchmark.mls.uzh.ch:9000",
            "access_key": "asdfasdf",
            "secret_key": "asdfasdfasdf",
            "secure": False,
        },
        benchmark="test",
    )
    ss._get_containers()
    assert ss.containers == ["test.0.1", "test.0.2", "test.1.0", "test2.0.1"]


def test__get_benchmarks(monkeypatch):
    monkeypatch.setattr(S3Storage, "connect", mock_test_connect)

    ss = S3Storage(
        auth_options={
            "endpoint": "http://omnibenchmark.mls.uzh.ch:9000",
            "access_key": "asdfasdf",
            "secret_key": "asdfasdfasdf",
            "secure": False,
        },
        benchmark="test",
    )
    ss._get_benchmarks()
    assert ss.benchmarks == ["test", "test2"]

    ss._get_benchmarks(update=False)
    assert ss.benchmarks == ["test", "test2"]


def test__create_benchmark(monkeypatch):
    monkeypatch.setattr(S3Storage, "connect", mock_test_connect)

    ss = S3Storage(
        auth_options={
            "endpoint": "http://omnibenchmark.mls.uzh.ch:9000",
            "access_key": "asdfasdf",
            "secret_key": "asdfasdfasdf",
            "secure": False,
        },
        benchmark="test",
    )
    ss._create_benchmark("test3")
    # assert ss.benchmarks == ['test', 'test2','test3']

    ss = S3Storage(
        auth_options={
            "endpoint": "http://omnibenchmark.mls.uzh.ch:9000",
            "access_key": "asdfasdf",
            "secret_key": "asdfasdfasdf",
            "secure": False,
        },
        benchmark="test",
    )
    ss._create_benchmark("test3", update=False)


def test__get_versions(monkeypatch):
    monkeypatch.setattr(S3Storage, "connect", mock_test_connect)

    ss = S3Storage(
        auth_options={
            "endpoint": "http://omnibenchmark.mls.uzh.ch:9000",
            "access_key": "asdfasdf",
            "secret_key": "asdfasdfasdf",
            "secure": False,
        },
        benchmark="test",
    )
    ss._get_versions()
    assert ss.versions == ["0.1", "0.2", "1.0"]

    ss._get_versions(update=False)
    assert ss.versions == ["0.1", "0.2", "1.0"]


def test__get_versions_public(monkeypatch):
    monkeypatch.setattr(requests, "get", mock_test_requests_get)

    ss = S3Storage(
        auth_options={
            "endpoint": "http://omnibenchmark.mls.uzh.ch:9000",
            "secure": False,
        },
        benchmark="test",
    )
    ss._get_versions()
    assert ss.versions == ["0.1", "0.2"]

    ss._get_versions(readonly=True)
    assert ss.versions == ["0.1", "0.2"]


def test__create_new_version(monkeypatch):
    monkeypatch.setattr(S3Storage, "connect", mock_test_connect)

    ss = S3Storage(
        auth_options={
            "endpoint": "http://omnibenchmark.mls.uzh.ch:9000",
            "access_key": "asdfasdf",
            "secret_key": "asdfasdfasdf",
            "secure": False,
        },
        benchmark="test",
    )
    ss.set_new_version()
    assert ss.major_version == 1
    assert ss.minor_version == 0
    assert ss.major_version_new == 1
    assert ss.minor_version_new == 1
    assert ss.versions == ["0.1", "0.2", "1.0"]

    # fails because no new version is created in mock test
    with pytest.raises(ValueError):
        ss._create_new_version()


def test__set_current_version(monkeypatch):
    monkeypatch.setattr(S3Storage, "connect", mock_test_connect)

    ss = S3Storage(
        auth_options={
            "endpoint": "http://omnibenchmark.mls.uzh.ch:9000",
            "access_key": "asdfasdf",
            "secret_key": "asdfasdfasdf",
            "secure": False,
        },
        benchmark="test",
    )
    ss.set_current_version()
    assert ss.major_version == 1
    assert ss.minor_version == 0


def test__set_current_version_public(monkeypatch):
    monkeypatch.setattr(requests, "get", mock_test_requests_get)

    ss = S3Storage(
        auth_options={
            "endpoint": "http://omnibenchmark.mls.uzh.ch:9000",
            "secure": False,
        },
        benchmark="test",
    )
    ss.set_current_version()
    assert ss.major_version == 0
    assert ss.minor_version == 2


def test__set_new_version(monkeypatch):
    monkeypatch.setattr(S3Storage, "connect", mock_test_connect)

    ss = S3Storage(
        auth_options={
            "endpoint": "http://omnibenchmark.mls.uzh.ch:9000",
            "access_key": "asdfasdf",
            "secret_key": "asdfasdfasdf",
            "secure": False,
        },
        benchmark="test",
    )
    ss.set_new_version()
    assert ss.major_version_new == 1
    assert ss.minor_version_new == 1


def test__get_objects(monkeypatch):
    monkeypatch.setattr(S3Storage, "connect", mock_test_connect)

    ss = S3Storage(
        auth_options={
            "endpoint": "http://omnibenchmark.mls.uzh.ch:9000",
            "access_key": "asdfasdf",
            "secret_key": "asdfasdfasdf",
            "secure": False,
        },
        benchmark="test",
    )
    ss.set_current_version()
    ss._get_objects()
    assert ss.files.keys() == {"file1.txt", "file2.txt"}
    assert ss.files["file1.txt"].keys() == {
        "hash",
        "last_modified",
        "size",
        "symlink_path",
    }
    assert ss.files["file2.txt"].keys() == {
        "hash",
        "last_modified",
        "size",
        "symlink_path",
    }
    assert ss.files["file1.txt"]["size"] == 1
    assert ss.files["file2.txt"]["size"] == 1
    assert ss.files["file1.txt"]["hash"] == "hash"
    assert ss.files["file2.txt"]["hash"] == "hash"
    assert ss.files["file1.txt"]["last_modified"] == "2021-07-08T14:33:00.000000"
    assert ss.files["file2.txt"]["last_modified"] == "2021-07-08T14:33:00.000000"


def test__get_objects_public(monkeypatch):
    monkeypatch.setattr(requests, "get", mock_test_requests_get)

    ss = S3Storage(
        auth_options={
            "endpoint": "http://omnibenchmark.mls.uzh.ch:9000",
            "secure": False,
        },
        benchmark="test",
    )
    with pytest.raises(ValueError):
        ss._get_objects()

    ss.set_current_version()
    ss._get_objects()
    assert ss.files.keys() == {"test.txt", "test2.txt"}
    assert ss.files["test.txt"].keys() == {
        "hash",
        "last_modified",
        "size",
        "x-object-meta-mtime",
        "symlink_path",
        "accesstime",
    }


def test_find_objects_to_copy(monkeypatch):
    monkeypatch.setattr(S3Storage, "connect", mock_test_connect)

    ss = S3Storage(
        auth_options={
            "endpoint": "http://omnibenchmark.mls.uzh.ch:9000",
            "access_key": "asdfasdf",
            "secret_key": "asdfasdfasdf",
            "secure": False,
        },
        benchmark="test",
    )
    ss.set_current_version()
    ss.find_objects_to_copy(tagging_type="all")
    assert ss.files.keys() == {"file1.txt", "file2.txt"}
    assert ss.files["file1.txt"].keys() == {
        "hash",
        "last_modified",
        "size",
        "symlink_path",
        "copy",
    }
    assert ss.files["file2.txt"].keys() == {
        "hash",
        "last_modified",
        "size",
        "symlink_path",
        "copy",
    }
    assert type(ss.files["file1.txt"]["copy"]) == bool
    assert type(ss.files["file2.txt"]["copy"]) == bool

    with pytest.raises(ValueError):
        ss.find_objects_to_copy(tagging_type="other")

    with pytest.raises(ValueError):
        ss.find_objects_to_copy(reference_time="2024")

    import datetime

    ss.find_objects_to_copy(reference_time=datetime.datetime.now())
    assert ss.files.keys() == {"file1.txt", "file2.txt"}
    assert ss.files["file1.txt"].keys() == {
        "hash",
        "last_modified",
        "size",
        "symlink_path",
        "copy",
    }
    assert ss.files["file2.txt"].keys() == {
        "hash",
        "last_modified",
        "size",
        "symlink_path",
        "copy",
    }
    assert ss.files["file1.txt"]["copy"] == True
    assert ss.files["file2.txt"]["copy"] == True


def test_copy_objects(monkeypatch):
    monkeypatch.setattr(S3Storage, "connect", mock_test_connect)

    ss = S3Storage(
        auth_options={
            "endpoint": "http://omnibenchmark.mls.uzh.ch:9000",
            "access_key": "asdfasdf",
            "secret_key": "asdfasdfasdf",
            "secure": False,
        },
        benchmark="test",
    )
    ss.set_new_version()
    ss.find_objects_to_copy()
    ss.copy_objects()
    assert ss.files.keys() == {"file1.txt", "file2.txt"}
    assert ss.files["file1.txt"].keys() == {
        "hash",
        "last_modified",
        "size",
        "symlink_path",
        "copy",
        "copied",
    }
    assert ss.files["file2.txt"].keys() == {
        "hash",
        "last_modified",
        "size",
        "symlink_path",
        "copy",
        "copied",
    }
    assert type(ss.files["file1.txt"]["copied"]) == bool
    assert type(ss.files["file2.txt"]["copied"]) == bool

    with pytest.raises(NotImplementedError):
        ss.copy_objects(type="symlink")
    with pytest.raises(ValueError):
        ss.copy_objects(type="other")


def test_create_new_version(monkeypatch):
    monkeypatch.setattr(S3Storage, "connect", mock_test_connect)

    ss = S3Storage(
        auth_options={
            "endpoint": "http://omnibenchmark.mls.uzh.ch:9000",
            "access_key": "asdfasdf",
            "secret_key": "asdfasdfasdf",
            "secure": False,
        },
        benchmark="test",
    )

    # fails because no new version is created in mock test
    with pytest.raises(ValueError):
        ss.create_new_version()
