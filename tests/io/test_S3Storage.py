import pytest
import re
from omni.io.S3Storage import S3Storage
import boto3
import json
import os
from pathlib import Path

if Path("tests/.minio_test_config.json").exists():
    with open("tests/.minio_test_config.json", "r") as file:
        auth_options = json.load(file)
elif (
    "S3_ENDPOINT_URL" in os.environ
    and "S3_ACCESS_KEY" in os.environ
    and "S3_SECRET_KEY" in os.environ
):
    auth_options = {}
    auth_options["endpoint"] = os.environ["S3_ENDPOINT_URL"]
    auth_options["access_key"] = os.environ["S3_ACCESS_KEY"]
    auth_options["secret_key"] = os.environ["S3_SECRET_KEY"]
    auth_options["secure"] = False
else:
    raise ValueError("No S3 credentials found")

auth_options_readonly = auth_options.copy()
del auth_options_readonly["access_key"]
del auth_options_readonly["secret_key"]


def cleanup_buckets():
    cleanupss = boto3.resource(
        "s3",
        endpoint_url=auth_options["endpoint"],
        aws_access_key_id=auth_options["access_key"],
        aws_secret_access_key=auth_options["secret_key"],
        aws_session_token=None,
        config=boto3.session.Config(signature_version="s3v4"),
        verify=False,
    )
    buckets = [i.name for i in cleanupss.buckets.all()]
    for bucket in buckets:
        if re.search(r"^test\.", bucket) or re.search(r"^test\d\.", bucket):
            cleanupss.Bucket(bucket).objects.all().delete()
            cleanupss.Bucket(bucket).delete()


@pytest.fixture(scope="session", autouse=True)
def cleanup(request):
    """Cleanup a testing directory once we are finished."""
    request.addfinalizer(cleanup_buckets)


def test_init():
    ss = S3Storage(auth_options=auth_options, benchmark="test")


def test__test_connect():
    ss = S3Storage(auth_options=auth_options, benchmark="test")
    assert ss.benchmark == "test"
    ss._test_connect()


def test_init_public():
    ss = S3Storage(auth_options=auth_options_readonly, benchmark="test")


def test__get_containers():
    ss = S3Storage(auth_options=auth_options, benchmark="test")
    ss._get_containers()
    assert [ct in ss.containers for ct in ["test.0.1"]] == [True]
    # assert ss.containers == ['test.0.1', 'test.0.2', 'test.1.0', 'test2.0.1']


def test__get_benchmarks():
    ss = S3Storage(auth_options=auth_options, benchmark="test")
    ss._get_benchmarks()
    assert "test" in ss.benchmarks

    ss._get_benchmarks(update=False)
    assert "test" in ss.benchmarks


def test__create_benchmark():
    cleanup_buckets()
    ss = S3Storage(auth_options=auth_options, benchmark="test")
    ss._create_benchmark("test2")
    # assert ss.benchmarks == ['test', 'test2','test3']

    cleanup_buckets()
    ss = S3Storage(auth_options=auth_options, benchmark="test")
    ss._create_benchmark("test3", update=False)


def test__get_versions():
    ss = S3Storage(auth_options=auth_options, benchmark="test")
    ss._get_versions()
    assert ss.versions == ["0.1"]

    ss._get_versions(update=False)
    assert ss.versions == ["0.1"]


def test__get_versions_public():
    ss = S3Storage(auth_options=auth_options_readonly, benchmark="test")
    ss._get_versions()
    assert ss.versions == ["0.1"]

    ss._get_versions(readonly=True)
    assert ss.versions == ["0.1"]


def test__create_new_version():
    ss = S3Storage(auth_options=auth_options, benchmark="test")
    ss.set_new_version()
    assert ss.major_version == 0
    assert ss.minor_version == 1
    assert ss.major_version_new == 0
    assert ss.minor_version_new == 2
    assert ss.versions == ["0.1"]


def test__set_current_version():
    ss = S3Storage(auth_options=auth_options, benchmark="test")
    ss.set_current_version()
    assert ss.major_version == 0
    assert ss.minor_version == 1


def test__set_current_version_public():
    ss = S3Storage(auth_options=auth_options_readonly, benchmark="test")
    ss.set_current_version()
    assert ss.major_version == 0
    assert ss.minor_version == 1


def test__set_new_version():
    ss = S3Storage(auth_options=auth_options, benchmark="test")
    ss.set_new_version()
    assert ss.major_version_new == 0
    assert ss.minor_version_new == 2


def test__get_objects():
    ss = S3Storage(auth_options=auth_options, benchmark="test")
    result = ss.client.Bucket(f"{ss.benchmark}.0.1").put_object(
        Key="file1.txt", Body=b""
    )
    result = ss.client.Bucket(f"{ss.benchmark}.0.1").put_object(
        Key="file2.txt", Body=b""
    )
    ss.set_current_version()
    ss._get_objects()
    print(ss.files)
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
    assert ss.files["file1.txt"]["size"] == 0
    assert ss.files["file2.txt"]["size"] == 0
    assert ss.files["file1.txt"]["hash"] == "d41d8cd98f00b204e9800998ecf8427e"
    assert ss.files["file2.txt"]["hash"] == "d41d8cd98f00b204e9800998ecf8427e"
    assert len(ss.files["file1.txt"]["last_modified"]) > 0
    assert len(ss.files["file2.txt"]["last_modified"]) > 0


def test__get_objects_public():
    ss = S3Storage(auth_options=auth_options_readonly, benchmark="test")
    # with pytest.raises(ValueError):
    #     ss._get_objects()

    ss.set_current_version()
    ss._get_objects()
    print(ss.files)
    assert ss.files.keys() == {"file1.txt", "file2.txt"}
    assert ss.files["file1.txt"].keys() == {
        "hash",
        "last_modified",
        "size",
        "x-object-meta-mtime",
        "symlink_path",
        "accesstime",
    }


def test_find_objects_to_copy():
    ss = S3Storage(auth_options=auth_options, benchmark="test")
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


def test_copy_objects():
    ss = S3Storage(auth_options=auth_options, benchmark="test")
    ss.set_new_version()
    ss.find_objects_to_copy()
    ss._create_new_version()
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


def test_create_new_version():
    ss = S3Storage(auth_options=auth_options, benchmark="test")

    ss.create_new_version()
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
