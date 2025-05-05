import zipfile
from pathlib import Path
from typing import List
import os

from omni_schema.datamodel import omni_schema

from omnibenchmark.benchmark import Benchmark
from omnibenchmark.io.code import clone_module
from omnibenchmark.io.files import download_files, list_files


def prepare_archive_code(benchmark: Benchmark) -> List[Path]:
    """
    Prepare the code to archive and return list of all filenames.

    Args:
        benchmark (str): The benchmark id.

    Returns:
        List[Path]: The filenames of all code to archive.
    """

    nodes = benchmark.get_nodes()
    repositories = set()
    for node in nodes:
        repositories.add((node.get_repository().url, node.get_repository().commit))

    repositories_dir = Path(".snakemake") / "repos"
    files = []
    while repositories:
        repo = repositories.pop()
        files += list(clone_module(repositories_dir, repo[0], repo[1]).iterdir())
    return files


def prepare_archive_software(benchmark: Benchmark) -> List[Path]:
    """
    Prepare the software to archive and return list of all filenames.

    Args:
        benchmark (str): The benchmark id.

    Returns:
        List[Path]: The filenames of all software to archive.
    """
    # decide on software type
    which_env = benchmark.get_benchmark_software_backend()
    files = []
    if which_env != omni_schema.SoftwareBackendEnum.host:
        if which_env == omni_schema.SoftwareBackendEnum.envmodules:
            files += prepare_archive_software_easyconfig(benchmark)
        elif which_env == omni_schema.SoftwareBackendEnum.apptainer:
            files += prepare_archive_software_apptainer(benchmark)
        elif which_env == omni_schema.SoftwareBackendEnum.conda:
            files += prepare_archive_software_conda(benchmark)
        else:
            raise NotImplementedError(
                f"Software backend {which_env} not implemented yet."
            )
    return files


def prepare_archive_software_easyconfig(benchmark: Benchmark) -> List[Path]:
    """
    Prepare the software easyconfig to archive and return list of all filenames.

    Args:
        benchmark (str): The benchmark id.

    Returns:
        List[Path]: The filenames of all software easyconfig to archive.
    """
    files = []
    softenvs = benchmark.get_benchmark_software_environments()
    for softenv in softenvs.values():
        envmodule_file = (benchmark.directory / Path(softenv.envmodule)).relative_to(
            Path(os.getcwd())
        )
        if envmodule_file.is_file():
            files.append(envmodule_file)
        else:
            raise FileNotFoundError(f"File {envmodule_file} not found.")
        easyconfig_file = (
            benchmark.directory / Path(softenv.easyconfig_file)
        ).relative_to(Path(os.getcwd()))
        if easyconfig_file.is_file():
            files.append(easyconfig_file)
        else:
            raise FileNotFoundError(f"File {easyconfig_file} not found.")
    return files


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
    # from omnibenchmark.software.conda_backend import pin_conda_envs
    # pin_conda_envs(benchmark.get_definition_file())
    files = []
    softenvs = benchmark.get_benchmark_software_environments()
    for softenv in softenvs.values():
        conda_file = (benchmark.directory / Path(softenv.conda)).relative_to(
            Path(os.getcwd())
        )
        if conda_file.is_file():
            files.append(conda_file)
        else:
            raise FileNotFoundError(f"File {conda_file} not found.")
    return files


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
    files = []
    softenvs = benchmark.get_benchmark_software_environments()
    for softenv in softenvs.values():
        apptainer_file = (benchmark.directory / Path(softenv.apptainer)).relative_to(
            Path(os.getcwd())
        )
        if apptainer_file.is_file():
            files.append(apptainer_file)
        else:
            raise FileNotFoundError(f"File {apptainer_file} not found.")
    return files


def prepare_archive_results(benchmark: Benchmark, local: bool = False) -> List[Path]:
    """
    Prepare the results to archive and return list of all filenames.

    Args:
        benchmark (str): The benchmark id.

    Returns:
        List[Path]: The filenames of all results to archive.
    """
    # get all results, check if exist locally, otherwise download

    objectnames, etags = list_files(
        benchmark=benchmark.get_definition_file(),
        type="all",
        stage=None,
        module=None,
        file_id=None,
        local=local,
    )
    if not local:
        download_files(
            benchmark=benchmark.get_definition_file(),
            type="all",
            stage=None,
            module=None,
            file_id=None,
            overwrite=True,
        )
    return objectnames


def archive_version(
    benchmark: Benchmark,
    outdir: Path = Path(),
    config: bool = True,
    code: bool = False,
    software: bool = False,
    results: bool = False,
    compression=zipfile.ZIP_STORED,
    compresslevel: int = None,
    dry_run: bool = False,
    local: bool = False,
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
    #### save easyconfig files
    ### conda
    #### save explicit packages list
    ### apptainer
    #### save .sif file
    if software:
        filenames += prepare_archive_software(benchmark)

    ## results (results files)
    ### check if results match remote, if not download
    if results:
        filenames += prepare_archive_results(benchmark, local)

    if dry_run:
        return filenames
    else:
        match compression:
            case zipfile.ZIP_BZIP2:
                file_extension = ".bz2"
            case zipfile.ZIP_LZMA:
                file_extension = ".xz"
            case _:
                file_extension = ".zip"
        # save all files to zip archive
        outfile = f"{benchmark.get_benchmark_name()}_{benchmark.get_converter().get_version()}{file_extension}"
        with zipfile.ZipFile(
            outdir / outfile, "w", compression=compression, compresslevel=compresslevel
        ) as archive:
            for filename in filenames:
                archive.write(filename, filename)

        return outdir / outfile


# archive_version(benchmark, outdir=Path(), config=True, code=True, software=False, results=False)
