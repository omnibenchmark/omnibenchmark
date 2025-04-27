import os
import random
import shutil
import string
import time
from pathlib import Path

import minio
import yaml
from linkml_runtime.dumpers import yaml_dumper
from testcontainers.minio import MinioContainer

from omni.benchmark import Benchmark


class MinIOSetup:
    def __init__(self, init: bool = True) -> None:
        self.do_init = init
        if self.do_init:
            self.minio = MinioContainer(
                image="minio/minio:RELEASE.2024-06-13T22-53-53Z"
            )
            self.minio.start()
            self.bucket_names = [
                "b"
                + "".join(random.choices(string.ascii_lowercase + string.digits, k=20))
                for _ in range(1000)
            ]

            time.sleep(2)

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

            self.auth_options_readonly["endpoint"] = (
                "http://"
                + self.auth_options_readonly["endpoint"]
                .replace("http://", "")
                .replace("https://", "")
            )

            self.bucket_name = testcontainer.bucket_names.pop()
            self.out_dir = ""

            os.environ["OB_STORAGE_S3_ACCESS_KEY"] = self.auth_options["access_key"]
            os.environ["OB_STORAGE_S3_SECRET_KEY"] = self.auth_options["secret_key"]

    def setup(self, in_dir: Path, out_dir: Path):
        if self.do_init:
            if not os.path.exists(out_dir):
                os.mkdir(out_dir)
            self.out_dir = out_dir / self.bucket_name
            if not os.path.exists(self.out_dir):
                os.mkdir(self.out_dir)

            filename_out = "Benchmark_" + self.bucket_name + ".yaml"
            with open(in_dir / "Benchmark_004.yaml", "r") as fh:
                yaml.safe_load(fh)
                benchmark_obj = Benchmark(Path(in_dir / "Benchmark_004.yaml"))

            benchmark_obj.converter.model.storage = self.auth_options["endpoint"]
            benchmark_obj.converter.model.storage_bucket_name = self.bucket_name
            yaml_dumper.dump(benchmark_obj.converter.model, self.out_dir / filename_out)
            if not os.path.exists(self.out_dir / "envs"):
                os.mkdir(self.out_dir / "envs")
            env_path = in_dir / "envs" / "python_vX_test.yaml"
            env_path_after = self.out_dir / "envs" / "python_vX_test.yaml"

            if not os.path.isfile(env_path_after):
                shutil.copyfile(env_path, env_path_after)

            self.benchmark_file = self.out_dir / filename_out

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
            if os.path.isdir(self.out_dir):
                shutil.rmtree(self.out_dir)
