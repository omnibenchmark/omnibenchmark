import re

from testcontainers.minio import MinioContainer

from omni.io.MinIOStorage import set_bucket_public_readonly


class MinIOSetup:
    def __init__(self, init: bool = True) -> None:
        self.do_init = init
        if self.do_init:
            self.minio = MinioContainer(
                image="minio/minio:RELEASE.2024-06-13T22-53-53Z"
            )
            self.minio.start()
            self.minio.get_client().make_bucket("benchmarks")
            set_bucket_public_readonly(self.minio.get_client(), "benchmarks")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if self.do_init:
            self.minio.stop()


class TmpMinIOStorage:
    def __init__(self, testcontainer) -> None:
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

        self.auth_options_readonly["endpoint"] = "http://" + self.auth_options_readonly[
            "endpoint"
        ].replace("http://", "").replace("https://", "")

        self.bucket_base = "test1"

    def cleanup_buckets(self):
        tmp_auth_options = self.auth_options.copy()
        tmp_auth_options["endpoint"] = (
            tmp_auth_options["endpoint"].replace("http://", "").replace("https://", "")
        )
        cleanupss = self.minio.get_client()
        buckets = [i.name for i in cleanupss.list_buckets()]
        for bucket in buckets:
            if re.search(rf"^{self.bucket_base}\.", bucket) or re.search(
                rf"^{self.bucket_base}\d\.", bucket
            ):
                objects = cleanupss.list_objects(bucket, recursive=True)
                for object in objects:
                    cleanupss.remove_object(bucket, object.object_name)
                cleanupss.remove_bucket(bucket)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.cleanup_buckets()
