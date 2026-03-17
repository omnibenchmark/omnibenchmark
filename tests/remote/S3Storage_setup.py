import os
import random
import shutil
import string
import time
import urllib.error
import urllib.request
from pathlib import Path

from testcontainers.core.container import DockerContainer

import yaml

from omnibenchmark.model import Benchmark
from omnibenchmark.remote.S3Storage import S3CompatibleStorage
from omnibenchmark.remote.RemoteStorage import StorageOptions

RUSTFS_IMAGE = "rustfs/rustfs:1.0.0-alpha.85"
RUSTFS_ACCESS_KEY = "rustfsadmin"
RUSTFS_SECRET_KEY = "rustfsadmin"

# How many buckets to pre-seed
BUCKETS_NUMBER = 100


class RustFSContainer:
    def __init__(
        self,
        image: str = RUSTFS_IMAGE,
        access_key: str = RUSTFS_ACCESS_KEY,
        secret_key: str = RUSTFS_SECRET_KEY,
    ) -> None:
        self.access_key = access_key
        self.secret_key = secret_key
        self._container = (
            DockerContainer(image)
            .with_exposed_ports(9000)
            .with_env("RUSTFS_ACCESS_KEY", self.access_key)
            .with_env("RUSTFS_SECRET_KEY", self.secret_key)
            .with_env("RUSTFS_VOLUMES", "/data")
        )

    def start(self):
        self._container.start()
        self._wait_for_ready()
        return self

    def stop(self):
        self._container.stop()

    def _wait_for_ready(self, timeout: int = 60):
        port = self._container.get_exposed_port(9000)
        url = f"http://localhost:{port}"
        deadline = time.time() + timeout
        while time.time() < deadline:
            try:
                urllib.request.urlopen(url, timeout=2)
                return
            except urllib.error.HTTPError:
                # Server is up but returned an HTTP error (e.g. 403) — that's fine
                return
            except Exception:
                time.sleep(1)
        raise TimeoutError(f"RustFS did not become ready within {timeout}s")

    def get_config(self):
        port = self._container.get_exposed_port(9000)
        return {"endpoint": f"localhost:{port}"}


class RustFSSetup:
    """Manages a RustFS container and pre-seeds buckets for testing."""

    def __init__(self) -> None:
        self.container = RustFSContainer()
        self.container.start()

        # Pre-seed buckets
        self.bucket_names = [
            "b" + "".join(random.choices(string.ascii_lowercase + string.digits, k=20))
            for _ in range(BUCKETS_NUMBER)
        ]

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.container.stop()


class TmpRustFSStorage:
    """Manages a temporary bucket for a single test."""

    def __init__(self, rustfs: RustFSSetup) -> None:
        self.container = rustfs.container
        self.endpoint = self.container.get_config()["endpoint"].replace(
            "localhost", "http://localhost"
        )
        self.access_key = self.container.access_key
        self.secret_key = self.container.secret_key

        # Use one of the pre-seeded buckets
        self.bucket_name = rustfs.bucket_names.pop()
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
        benchmark_obj = Benchmark.from_yaml(Path(in_dir / benchmark_file))
        # Update the storage configuration in the benchmark model
        from omnibenchmark.model import Storage, StorageAPIEnum

        benchmark_obj.storage = Storage(
            api=StorageAPIEnum.s3, endpoint=self.endpoint, bucket_name=self.bucket_name
        )
        benchmark_file = str(self.out_dir / f"Benchmark_{self.bucket_name}.yaml")
        self.benchmark_file = benchmark_file
        # Use Pydantic's model dump with mode='json' to properly serialize enums
        # Exclude deprecated storage fields to avoid "mixed format" validation errors
        with open(benchmark_file, "w") as f:
            yaml.dump(
                benchmark_obj.model_dump(
                    mode="json", exclude={"storage_api", "storage_bucket_name"}
                ),
                f,
            )

        self.storage_options = StorageOptions(out_dir="out")

        # Copy environment file
        env_path = in_dir / "envs" / "python_vX_test.yaml"
        env_path_after = self.out_dir / "envs" / "python_vX_test.yaml"
        os.makedirs(self.out_dir / "envs", exist_ok=True)
        if not os.path.isfile(env_path_after):
            shutil.copyfile(env_path, env_path_after)

    def get_storage_client(self):
        return S3CompatibleStorage(
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
