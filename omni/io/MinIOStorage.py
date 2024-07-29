"""MinIO class for remote storage."""

import datetime
import io
import json
import logging
import os
import re
from typing import Dict
from urllib.parse import urlparse

import dateutil.parser
import minio
import minio.commonconfig
import minio.deleteobjects
import minio.retention
import requests
from packaging.version import Version

from omni.io.RemoteStorage import RemoteStorage, is_valid_version
from omni.io.S3config import S3_DEFAULT_STORAGE_OPTIONS, bucket_readonly_policy

logging.basicConfig(level=logging.ERROR)
logging.getLogger("requests").setLevel(logging.DEBUG)
logging.getLogger("minio").setLevel(logging.DEBUG)
logger = logging.getLogger(__name__)


def get_s3_object_versions_and_tags(client, benchmark, readonly=False):
    if not client.bucket_exists(benchmark):
        raise ValueError(f"Benchmark {benchmark} does not exist.")

    object_ls = list(
        client.list_objects(
            benchmark, include_version=True, include_user_meta=True, recursive=True
        )
    )
    object_names = [x.object_name for x in object_ls]
    version_ids = [x.version_id for x in object_ls]
    sizes = [x.size for x in object_ls]
    last_modified = [x.last_modified for x in object_ls]
    uq_names = list(set(object_names))

    di = {}
    for uq in uq_names:
        di[uq] = {}
        # tv = [v for v,n in zip(version_ids,object_names) if n == uq]
        tmpinds = [i for i, n in enumerate(object_names) if n == uq]
        tv = [version_ids[i] for i in tmpinds]
        ts = [sizes[i] for i in tmpinds]
        tlm = [last_modified[i] for i in tmpinds]
        tdl = [object_ls[i].is_delete_marker for i in tmpinds]
        etl = [object_ls[i].etag for i in tmpinds]
        for j, v in enumerate(tv):
            tags_filt = {}
            if not readonly:
                if not tdl[j]:
                    tags = client.get_object_tags(benchmark, uq, version_id=v)
                    if not tags is None:
                        for k, w in tags.items():
                            if is_valid_version(k):
                                tags_filt[k] = w
            di[uq][v] = {}
            di[uq][v]["tags"] = tags_filt
            tmpinds2 = [i for i, tvers in enumerate(tv) if tvers == v]
            di[uq][v]["size"] = ts[tmpinds2[0]]
            di[uq][v]["last_modified"] = tlm[tmpinds2[0]]
            di[uq][v]["is_delete_marker"] = tdl[tmpinds2[0]]
            di[uq][v]["etag"] = etl[tmpinds2[0]]
    return di


def get_s3_objects_to_tag(
    objdic, tracked_dirs=S3_DEFAULT_STORAGE_OPTIONS["tracked_directories"]
):
    object_names = sorted(list(objdic.keys()))
    newest_versions = []
    for uq in object_names:
        subdatetimels = []
        for v in objdic[uq].keys():
            subdatetimels.append(objdic[uq][v]["last_modified"])
        maxind = [k for k, sd in enumerate(subdatetimels) if sd == max(subdatetimels)][
            0
        ]
        newest_versions.append(list(objdic[uq].keys())[maxind])

    newest_version_exists = []
    for uq, v in zip(object_names, newest_versions):
        newest_version_exists.append(not objdic[uq][v]["is_delete_marker"])

    is_tracked_basedirls = []
    for n in object_names:
        # get root directory
        while not os.path.split(n)[0] == "":
            n = os.path.split(n)[0]
        if n in tracked_dirs:
            is_tracked_basedirls.append(True)
        else:
            is_tracked_basedirls.append(False)

    do_tag = []
    for nve, itb in zip(newest_version_exists, is_tracked_basedirls):
        if nve and itb:
            do_tag.append(True)
        else:
            do_tag.append(False)
    object_names_to_tag = [n for n, dt in zip(object_names, do_tag) if dt]
    versionid_of_objects_to_tag = [
        newest_versions[i] for i, dt in enumerate(do_tag) if dt
    ]
    return object_names_to_tag, versionid_of_objects_to_tag


def get_single_s3version_from_bmversion(di, object_name, query_version):
    return list(
        filter(
            lambda v: query_version in di[object_name][v]["tags"].keys(),
            di[object_name].keys(),
        )
    )


def get_s3version_from_bmversion(di, query_version):
    vv_ls = []
    for uq in di.keys():
        tmpv = get_single_s3version_from_bmversion(di, uq, query_version)
        if len(tmpv) == 0:
            continue
        elif len(tmpv) > 1:
            raise ValueError("Multiple versions found")
        else:
            vv_ls.append(
                [
                    uq,
                    tmpv[0],
                    di[uq][tmpv[0]]["last_modified"],
                    di[uq][tmpv[0]]["size"],
                    di[uq][tmpv[0]]["etag"],
                ]
            )
    return vv_ls


def prepare_csv_s3version_from_bmversion(vv_ls):
    outstr = "name,version_id,last_modified,size,etag\n"
    for vv in vv_ls:
        outstr += f"{vv[0]},{vv[1]},{vv[2]},{vv[3]},{vv[4]}\n"
    return outstr


def get_meta_mtime(preauthurl, containername, objectname):
    """
    Retrieves the metadata modification time and access time of a file from a MinIO storage.

    Args:
        preauthurl (str): The pre-authenticated URL of the MinIO storage.
        containername (str): The name of the container where the file is stored.
        objectname (str): The name of the file.

    Returns:
        tuple: A tuple containing the file's metadata modification time and access time.

    Raises:
        requests.HTTPError: If the HTTP request to retrieve the file fails.
    """
    urlfile = f"{preauthurl}/{containername}/{objectname}"
    response = requests.get(urlfile)
    if response.ok:
        response_headers = response.headers
        if "X-Object-Meta-Mtime" in response_headers.keys():
            file_time = datetime.datetime.fromtimestamp(
                float(response_headers["X-Object-Meta-Mtime"]), datetime.timezone.utc
            )
        elif "X-Object-Meta-Last-Modified" in response_headers.keys():
            file_time = dateutil.parser.parse(
                response_headers["X-Object-Meta-Last-Modified"]
            )
        else:
            file_time = dateutil.parser.parse(response_headers["Last-Modified"])
        access_time = dateutil.parser.parse(response_headers["Date"])
        return file_time, access_time
    else:
        response.raise_for_status()


def set_bucket_public_readonly(client, bucket_name):
    policy = bucket_readonly_policy(bucket_name)
    client.set_bucket_policy(bucket_name, json.dumps(policy))


def set_bucket_lifecycle_config(client, bucket_name, noncurrent_days=1):
    lifecycle_config = minio.lifecycleconfig.LifecycleConfig(
        [
            minio.lifecycleconfig.Rule(
                minio.commonconfig.ENABLED,
                rule_filter=minio.lifecycleconfig.Filter(prefix="*"),
                rule_id="rule1",
                noncurrent_version_expiration=minio.lifecycleconfig.NoncurrentVersionExpiration(
                    noncurrent_days=noncurrent_days
                ),
            ),
        ],
    )
    client.set_bucket_lifecycle(bucket_name, lifecycle_config)


class MinIOStorage(RemoteStorage):
    def __init__(
        self,
        auth_options: Dict,
        benchmark: str,
        storage_options: Dict = S3_DEFAULT_STORAGE_OPTIONS,
    ):
        super().__init__(auth_options, benchmark)
        assert "endpoint" in self.auth_options.keys()
        if "access_key" in self.auth_options.keys():
            self.client = self.connect()
            self.roclient = self.connect(readonly=True)
            self._test_connect()

            if not self.client.bucket_exists(benchmark):
                logger.warning(
                    f"Benchmark {benchmark} does not exist, creating new benchmark."
                )
                self._create_benchmark(benchmark)

            self._get_versions()
        else:
            self.roclient = self.connect(readonly=True)
            self._get_versions()
            # read-only mode
            pass

    def connect(self, readonly=False):
        """
        Connects to the MinIO storage.

        Returns:
        - A MinIO client object.
        """
        if readonly:
            assert "endpoint" in self.auth_options.keys()
            tmp_auth_options = self.auth_options.copy()
            if "access_key" in tmp_auth_options.keys():
                del tmp_auth_options["access_key"]
            if "secret_key" in tmp_auth_options.keys():
                del tmp_auth_options["secret_key"]
        else:
            assert (
                "endpoint" in self.auth_options.keys()
                and "access_key" in self.auth_options.keys()
                and "secret_key" in self.auth_options.keys()
            )
            tmp_auth_options = self.auth_options.copy()
        try:
            return minio.Minio(**tmp_auth_options)
        except Exception as e:
            url = urlparse(tmp_auth_options["endpoint"])
            tmp_auth_options["endpoint"] = url.netloc
            return minio.Minio(**tmp_auth_options)

    def _test_connect(self) -> None:
        try:
            _ = self.client.list_buckets()
        except Exception as e:
            raise e

    def _create_benchmark(self, benchmark):
        if self.client.bucket_exists(benchmark):
            raise ValueError("Benchmark already exists")
        # create new version
        self.client.make_bucket(bucket_name=benchmark, object_lock=True)
        set_bucket_public_readonly(self.client, benchmark)
        set_bucket_lifecycle_config(self.client, benchmark, noncurrent_days=1)
        _ = self.client.put_object(
            self.benchmark, "versions/0.0.csv", io.BytesIO(b""), 0
        )
        if not self.client.bucket_exists(benchmark):
            raise Exception(f"Bucket creation for benchmark {benchmark} failed")

    def _get_versions(self):
        versionobjects = list(
            self.roclient.list_objects(
                self.benchmark, prefix="versions", recursive=True
            )
        )
        allversions = [
            os.path.basename(v.object_name).replace(".csv", "") for v in versionobjects
        ]
        versions = list()
        for version in allversions:
            if re.match(f"(\\d+).(\\d+)", version):
                versions.append(Version(version))
        self.versions = versions

    def create_new_version(self):
        if self.version is None:
            raise ValueError(
                "No version provided, set version first with method 'set_version'"
            )

        if not "access_key" in self.auth_options.keys():
            raise ValueError("Read-only mode, cannot create new version")

        # update
        self._get_versions()
        if self.version in self.versions:
            raise ValueError("Version already exists")

        # get all objects
        objdic = get_s3_object_versions_and_tags(self.client, self.benchmark)

        # get objects to tag
        object_names_to_tag, versionid_of_objects_to_tag = get_s3_objects_to_tag(
            objdic, tracked_dirs=self.storage_options["tracked_directories"]
        )

        # Tag all objects with current version
        tags = minio.datatypes.Tags.new_object_tags()
        tags[str(self.version)] = "1"
        for n, v in zip(object_names_to_tag, versionid_of_objects_to_tag):
            self.client.set_object_tags(self.benchmark, n, tags, version_id=v)

        # set retention policy of objects
        retention_config = minio.retention.Retention(
            minio.commonconfig.GOVERNANCE,
            datetime.datetime.utcnow() + datetime.timedelta(weeks=1000),
        )
        for n, v in zip(object_names_to_tag, versionid_of_objects_to_tag):
            self.client.set_object_retention(
                self.benchmark, n, config=retention_config, version_id=v
            )

        # refresh object list
        objdic = get_s3_object_versions_and_tags(self.client, self.benchmark)

        # write versions/{{version}}.csv
        vv_ls = get_s3version_from_bmversion(objdic, str(self.version))
        vv_str = prepare_csv_s3version_from_bmversion(vv_ls)
        version_filename = f"versions/{self.version}.csv"
        _ = self.client.put_object(
            self.benchmark, version_filename, io.BytesIO(vv_str.encode()), len(vv_str)
        )

        # tag and set retention policy of version file
        self.client.set_object_tags(self.benchmark, version_filename, tags)
        self.client.set_object_retention(
            self.benchmark, version_filename, config=retention_config
        )

        # check if version exists
        self._get_versions()
        if not self.version in self.versions:
            raise ValueError("Version creation failed")

    def _get_objects(self, readonly=False):
        if self.version is None:
            raise ValueError("No version provided")

        if self.version in self.versions:
            # read overview file
            response = self.roclient.get_object(
                self.benchmark, f"versions/{self.version}.csv"
            )
            objls = response.data.decode("utf-8")
            objls = objls.split("\n")
            objls = [obj for obj in objls if obj]
            objdict = {}
            header = objls[0].split(",")
            assert header[0] == "name"
            for obj in objls[1:]:
                tmpsplit = obj.split(",")
                objdict[tmpsplit[0]] = {}
                for i, head in enumerate(header[1:]):
                    objdict[tmpsplit[0]][head] = tmpsplit[i + 1]
        else:
            # get all objects
            objdic = get_s3_object_versions_and_tags(
                self.roclient, self.benchmark, readonly=True
            )

            # get objects to tag
            object_names_to_tag, versionid_of_objects_to_tag = get_s3_objects_to_tag(
                objdic, tracked_dirs=self.storage_options["tracked_directories"]
            )
            objdict = {}
            for obj, vt in zip(object_names_to_tag, versionid_of_objects_to_tag):
                objdict[obj] = {}
                objdict[obj]["version_id"] = vt
                objdict[obj]["etag"] = objdic[obj][vt]["etag"]
                objdict[obj]["size"] = objdic[obj][vt]["size"]
                objdict[obj]["last_modified"] = objdic[obj][vt]["last_modified"]

        self.files = objdict

    def archive_version(self, version):
        NotImplementedError
        # self._update_overview(cleanup=True)

    def delete_version(self, version):
        NotImplementedError
        # self._update_overview(cleanup=True)


RemoteStorage.register(MinIOStorage)
