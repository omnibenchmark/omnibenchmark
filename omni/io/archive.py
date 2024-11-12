import zipfile
from pathlib import Path
from typing import List

from omni.benchmark import Benchmark


def prepare_archive_code(benchmark: Benchmark) -> List[Path]:
    """
    Prepare the code to archive and return list of all filenames.

    Args:
        benchmark (str): The benchmark id.

    Returns:
        List[Path]: The filenames of all code to archive.
    """
    # get all code repos (check local cache)
    # return all files
    return []


def prepare_archive_software(benchmark: Benchmark) -> List[Path]:
    """
    Prepare the software to archive and return list of all filenames.

    Args:
        benchmark (str): The benchmark id.

    Returns:
        List[Path]: The filenames of all software to archive.
    """
    # decide on software type
    # call prepare_archive_software_<type>
    return []


def prepare_archive_software_easyconfig(benchmark: Benchmark) -> List[Path]:
    """
    Prepare the software easyconfig to archive and return list of all filenames.

    Args:
        benchmark (str): The benchmark id.

    Returns:
        List[Path]: The filenames of all software easyconfig to archive.
    """
    # prepare easyconfigs software, get binaries, etc.
    # how?
    # return all files
    return []


def prepare_archive_software_envmodules(benchmark: Benchmark) -> List[Path]:
    """
    Prepare the envmodules to archive and return list of all filenames.

    Args:
        benchmark (str): The benchmark id.

    Returns:
        List[Path]: The filenames of all envmodules to archive.
    """
    # prepare envmodules software, get binaries, etc.
    # how?
    # return all files
    return []


def prepare_archive_software_conda(benchmark: Benchmark) -> List[Path]:
    """
    Prepare the conda to archive and return list of all filenames.

    Args:
        benchmark (str): The benchmark id.

    Returns:
        List[Path]: The filenames of all conda to archive.
    """
    # prepare conda, as in snakemake archive (https://github.com/snakemake/snakemake/blob/76d53290a003891c5ee41f81e8eb4821c406255d/snakemake/deployment/conda.py#L316) maybe?
    # return all files
    return []


def prepare_archive_software_apptainer(benchmark: Benchmark) -> List[Path]:
    """
    Prepare the software apptainer to archive and return list of all filenames.

    Args:
        benchmark (str): The benchmark id.

    Returns:
        List[Path]: The filenames of all software apptainer to archive.
    """
    # prepare apptainer, check if .sif file exists
    # return all files
    return []


def prepare_archive_results(benchmark: Benchmark) -> List[Path]:
    """
    Prepare the results to archive and return list of all filenames.

    Args:
        benchmark (str): The benchmark id.

    Returns:
        List[Path]: The filenames of all results to archive.
    """
    # get all results, check if exist locally, otherwise download
    return []


def archive_version(
    benchmark: Benchmark,
    outdir: Path = Path(),
    config: bool = True,
    code: bool = False,
    software: bool = False,
    results: bool = False,
):
    # retrieve all filenames to save
    filenames = []

    ## config (benchmark.yaml)
    if config:
        filenames.append(benchmark.get_definition_file())

    ## code (code files)
    ### check local cache of GH repos
    ### download repos if not in cache
    ### iterate over all files in repos
    if code:
        filenames += prepare_archive_code(benchmark)

    ## software (software files)
    ### easyconfig
    #### ?
    ### envmodule
    #### ?
    ### conda
    #### reuse snakemake archive (https://github.com/snakemake/snakemake/blob/76d53290a003891c5ee41f81e8eb4821c406255d/snakemake/deployment/conda.py#L316)
    ### apptainer
    #### save .sif file
    if software:
        filenames += prepare_archive_software(benchmark)

    ## results (results files)
    ### check if results match remote, if not download
    if results:
        filenames += prepare_archive_results(benchmark)

    # save all files to zip archive
    with zipfile.ZipFile(outdir / f"{benchmark.id}.zip", "w") as archive:
        for filename in filenames:
            archive.write(filename, filename)

    return outdir / f"{benchmark.id}.zip"
