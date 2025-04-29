from pathlib import Path


def get_benchmark_data_path() -> Path:
    return Path(__file__).resolve().parent.parent / "data"


data = get_benchmark_data_path()
