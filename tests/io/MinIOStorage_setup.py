import os
import random
import shutil
import string
from pathlib import Path

from testcontainers.minio import MinioContainer

# TODO(ben): deprecate in favor of pydantic model serializer
from linkml_runtime.dumpers import yaml_dumper

from omnibenchmark.benchmark import Benchmark
from omnibenchmark.io.MinIOStorage import MinIOStorage
from omnibenchmark.io.RemoteStorage import StorageOptions

MINIO_IMAGE = "minio/minio:RELEASE.2024-06-13T22-53-53Z"

# How many buckets to pre-seed
BUCKETS_NUMBER = 100


class MinIOSetup:
    """Manages a MinIO container and pre-seeds buckets for testing."""

    def __init__(self) -> None:
        self.minio = MinioContainer(image=MINIO_IMAGE)
        self.minio.start()

        # Pre-seed buckets
        self.bucket_names = [
            "b" + "".join(random.choices(string.ascii_lowercase + string.digits, k=20))
            for _ in range(BUCKETS_NUMBER)
        ]

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.minio.stop()


class TmpMinIOStorage:
    """Manages a temporary bucket for a single test."""

    def __init__(self, testcontainer: MinIOSetup) -> None:
        self.minio = testcontainer.minio
        self.endpoint = self.minio.get_config()["endpoint"].replace(
            "localhost", "http://localhost"
        )
        self.access_key = self.minio.access_key
        self.secret_key = self.minio.secret_key

        # Use one of the pre-seeded buckets
        self.bucket_name = testcontainer.bucket_names.pop()
        self.out_dir = ""

        # Set environment variables
        os.environ["OB_STORAGE_S3_ACCESS_KEY"] = self.access_key
        os.environ["OB_STORAGE_S3_SECRET_KEY"] = self.secret_key

    def setup(
        self, in_dir: Path, out_dir: Path, benchmark_file: str = "Benchmark_004.yaml"
    ):
        """Prepares the benchmark file and environment for testing."""
        self.out_dir = out_dir / self.bucket_name
        os.makedirs(self.out_dir, exist_ok=True)

        # Prepare benchmark file by injecting bucket name and endpoint
        benchmark_obj = Benchmark(Path(in_dir / benchmark_file))
        benchmark_obj.converter.model.storage = self.endpoint
        benchmark_obj.converter.model.storage_bucket_name = self.bucket_name
        benchmark_file = str(self.out_dir / f"Benchmark_{self.bucket_name}.yaml")
        self.benchmark_file = benchmark_file
        yaml_dumper.dump(benchmark_obj.converter.model, benchmark_file)

        self.storage_options = StorageOptions(out_dir="out")

        # Copy environment file
        env_path = in_dir / "envs" / "python_vX_test.yaml"
        env_path_after = self.out_dir / "envs" / "python_vX_test.yaml"
        os.makedirs(self.out_dir / "envs", exist_ok=True)
        if not os.path.isfile(env_path_after):
            shutil.copyfile(env_path, env_path_after)

    def get_storage_client(self):
        return MinIOStorage(
            auth_options=self.auth_options,
            storage_options=self.storage_options,
            benchmark=self.bucket_name,
        )

    @property
    def benchmark(self):
        return self.bucket_name

    @property
    def auth_options(self):
        return {
            "endpoint": self.endpoint,
            "access_key": self.access_key,
            "secret_key": self.secret_key,
            "secure": False,
        }

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        pass
