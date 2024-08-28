"""Sync benchmark definition file"""

import json
import os
import requests

from .config import bench_dir

base = "https://github.com/"
omni_essentials = "omnibenchmark/omni_essentials/-/raw/master/benchmarks/"
bench_definition = "benchmark_definition.yaml"


def get_bench_definition(bench_name: str, version: str, force: bool = False):
    bench_yaml = f"{bench_name}/{version}/{bench_definition}"
    local_bench_yaml = bench_dir + bench_yaml
    remote_bench_yaml = base + omni_essentials + f"{bench_yaml}?inline=false"
    if force or not os.path.isfile(local_bench_yaml):
        r = requests.get(remote_bench_yaml)
        data = r.json()
        with open(local_bench_yaml, "w") as f:
            json.dump(data, f)
