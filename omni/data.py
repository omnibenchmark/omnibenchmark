"""Functions to manage dataset handling"""


from .sync import get_bench_definition


def list_datafiles(bench_name: str, version: str, stage: str):
    """List all available datafiles for a certain benchmark, version and stage"""
    bench_yaml = get_bench_definition(bench_name, version, stage)


def get_datafile(url: str):
    """Download a specific datafile based on its url"""
    # Do we download/store/cache locally here?
    pass


def get_test_datafiles(bench_name: str, version: str, stage: str):
    """Download a test datafilesfor a certain benchmark version and stage"""
    pass
