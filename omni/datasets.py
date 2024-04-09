"""Functions to manage dataset handling"""


from .sync import get_bench_definition


def list_datasets(bench_name: str, version: str, stage: str):
    """List all available datasets for a certain benchmark, version and stage"""
    bench_yaml = get_bench_definition(bench_name, version, stage)


def mount_dataset(url: str):
    """Mount a specific dataset based on its url"""
    # Do we download/store/cache locally here?
    pass


def get_test_dataset(bench_name: str, version: str, stage: str):
    """Mount a test dataset for a certain benchmark version and stage"""
    pass