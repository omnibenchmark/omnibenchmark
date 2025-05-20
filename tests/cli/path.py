from pathlib import Path


def get_benchmark_data_path() -> Path:
    return Path(__file__).resolve().parent.parent / "data"


def get_benchmark_envs_path() -> Path:
    return Path(__file__).resolve().parent.parent / "envs"


data = get_benchmark_data_path()
