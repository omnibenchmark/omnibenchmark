"""Functions to manage dataset handling"""


from .sync import get_bench_definition


def list_output_files(bench_name: str, version: str, stage: str):
    """List all available files for a certain benchmark, version and stage"""
    bench_yaml = get_bench_definition(bench_name, version, stage)


def get_file(url: str):
    """Download a specific file based on its url"""
    # Do we download/store/cache locally here?
    pass


def get_test_inputs(bench_name: str, version: str, stage: str):
    """Download a test input for a certain benchmark version and stage"""
    pass
