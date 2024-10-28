import os
import re
import shutil
import tempfile
import time
from pathlib import Path
from typing import Union

import minio
import yaml
from linkml_runtime.dumpers import yaml_dumper
from testcontainers.minio import MinioContainer

from omni.benchmark import Benchmark
from omni.io.MinIOStorage import set_bucket_public_readonly


class MinIOSetup:
    def __init__(self, init: bool = True) -> None:
        self.do_init = init
        if self.do_init:
            self.minio = MinioContainer(
                image="minio/minio:RELEASE.2024-06-13T22-53-53Z"
            )
            self.minio.start()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if self.do_init:
            self.minio.stop()


class TmpMinIOStorage:
    def __init__(self, testcontainer) -> None:
        self.do_init = testcontainer.do_init
        if self.do_init:
            self.minio = testcontainer.minio
            self.auth_options = {}
            self.auth_options["endpoint"] = self.minio.get_config()["endpoint"].replace(
                "localhost", "http://localhost"
            )
            self.auth_options["access_key"] = self.minio.access_key
            self.auth_options["secret_key"] = self.minio.secret_key
            self.auth_options["secure"] = False

            self.auth_options_readonly = self.auth_options.copy()
            del self.auth_options_readonly["access_key"]
            del self.auth_options_readonly["secret_key"]

            self.auth_options_readonly[
                "endpoint"
            ] = "http://" + self.auth_options_readonly["endpoint"].replace(
                "http://", ""
            ).replace(
                "https://", ""
            )

            self.bucket_base = "test1"

    def create_yaml(self, in_dir, out_dir):
        if self.do_init:
            with open(in_dir / "Benchmark_004.yaml", "r") as fh:
                yaml.safe_load(fh)
                benchmark_obj = Benchmark(Path(in_dir / "Benchmark_004.yaml"))

            benchmark_obj.converter.model.storage = self.auth_options["endpoint"]
            yaml_dumper.dump(
                benchmark_obj.converter.model, out_dir / "Benchmark_004.yaml"
            )

            if not os.path.exists(out_dir / "envs"):
                os.mkdir(out_dir / "envs")
            env_path = in_dir / "envs" / "python_vX_test.yaml"
            env_path_after = out_dir / "envs" / "python_vX_test.yaml"

            shutil.copyfile(env_path, env_path_after)

            return out_dir / "Benchmark_004.yaml"

    def cleanup_buckets(self):
        if self.do_init:
            tmp_auth_options = self.auth_options.copy()
            tmp_auth_options["endpoint"] = (
                tmp_auth_options["endpoint"]
                .replace("http://", "")
                .replace("https://", "")
            )
            cleanupss = self.minio.get_client()
            # cleanupss = ss.client
            buckets = [i.name for i in cleanupss.list_buckets()]
            for bucket in buckets:
                objects = cleanupss.list_objects(
                    bucket, recursive=True, include_version=True
                )
                delobjectls = []
                for object in objects:
                    delobjectls.append(
                        minio.deleteobjects.DeleteObject(
                            name=object.object_name, version_id=object.version_id
                        )
                    )
                errors = cleanupss.remove_objects(
                    bucket, delobjectls, bypass_governance_mode=True
                )
                for error in errors:
                    print("error occurred when deleting object", error)
                cleanupss.remove_bucket(bucket)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if self.do_init:
            self.cleanup_buckets()


def create_remote_test(minio_testcontainer, in_dir, out_dir: Union[Path, None] = None):
    if out_dir is None:
        out_dir = Path(tempfile.gettempdir()) / "ob_test_benchmark004"
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    # need to create new benchmark yaml file with correct endpoint from test container
    with TmpMinIOStorage(minio_testcontainer) as tmp:
        time.sleep(2)
        os.environ["OB_STORAGE_S3_ACCESS_KEY"] = tmp.auth_options["access_key"]
        os.environ["OB_STORAGE_S3_SECRET_KEY"] = tmp.auth_options["secret_key"]
        return tmp.create_yaml(in_dir, out_dir)
