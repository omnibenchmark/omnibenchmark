import logging
import os
import subprocess
from pathlib import Path
from typing import List, Optional

import configparser

from snakemake.script import Snakemake

from omnibenchmark.benchmark.params import Params


def mock_execution(inputs: List[str], output: str, snakemake: Snakemake):
    print("Processed", inputs, "to", output, "using threads", snakemake.threads)
    print("  bench_iteration is", snakemake.bench_iteration)
    print("  resources are", snakemake.resources)
    print("  wildcards are", snakemake.wildcards)
    print("  rule is", snakemake.rule)
    print("  scriptdir is", snakemake.scriptdir)
    print("  params are", snakemake.params)


def execution(
    module_dir: Path,
    module_name: str,
    output_dir: Path,
    dataset: str,
    inputs_map: dict[str, str | List[str]],
    parameters: Params,
    keep_module_logs: bool,
    timeout: int,
) -> int:
    config_parser = _read_config(module_dir)
    if config_parser is None:
        logging.error(f"ERROR: Failed to load config.cfg for {module_name}.")
        return -1  # return exit_code -1

    executable = config_parser["DEFAULT"]["SCRIPT"]
    executable_path = module_dir / executable
    if not os.path.exists(executable_path):
        logging.error(f"ERROR: {module_name} executable does not exist.")
        return -1  # return exit_code -1

    # Constructing the command
    command = _create_command(executable_path)
    if command is None:
        return -1  # return exit_code -1

    # Add output directory and dataset name to command
    command.extend(["--output_dir", output_dir.as_posix(), "--name", dataset])

    # Adding input files with their respective keys
    if inputs_map:
        for k, v in inputs_map.items():
            if isinstance(v, str):
                command.extend([f"--{k}", Path(v).as_posix()])
            elif isinstance(v, list) and all(isinstance(item, str) for item in v):
                command.extend([f"--{k}"] + [Path(item).as_posix() for item in v])

    # Adding extra parameters
    if parameters:
        command.extend(parameters.to_cli_args())

    # Prepare stdout and stderr files in the output directory
    stdout_file = output_dir / "stdout.log"
    stderr_file = output_dir / "stderr.log"

    logging.info(f"Seting timeout to {timeout} seconds")

    try:
        with open(stdout_file, "w") as stdout_f, open(stderr_file, "w") as stderr_f:
            result = subprocess.run(
                command,
                check=True,
                stdout=stdout_f,
                stderr=stderr_f,
                text=True,
                timeout=timeout,
            )
            exit_code = result.returncode
            logging.info(f"{executable} ran successfully with return code {exit_code}.")

    except subprocess.CalledProcessError as e:
        exit_code = e.returncode
        logging.error(
            f"ERROR: Executing {executable} failed with exit code {e.returncode}. See {stderr_file} for details."
        )

    except subprocess.TimeoutExpired as e:
        logging.error(
            f"ERROR: Executing {executable} failed after timing out in {timeout}s"
            f"See {stderr_file} for details."
        )
        raise e

    finally:
        # Cleanup empty log files / always cleanup stdout_file if keep_module_logs is False
        if stdout_file.exists() and (
            os.path.getsize(stdout_file) == 0 or not keep_module_logs
        ):
            stdout_file.unlink()
        if stderr_file.exists() and os.path.getsize(stderr_file) == 0:
            stderr_file.unlink()

    if exit_code is None:
        logging.error(f"ERROR: Failed to execute {executable}: Unknown error")
        exit_code = -1

    return exit_code


def _read_config(module_dir: Path) -> Optional[configparser.ConfigParser]:
    config_path = module_dir / "config.cfg"
    if not os.path.exists(config_path):
        return None

    parser = configparser.ConfigParser()

    # Read the configuration file
    parser.read(config_path)

    return parser


def _create_command(executable_path: Path) -> Optional[list[str]]:
    _, extension = os.path.splitext(executable_path)

    if extension == ".py":
        # Execute Python script
        return ["python3", executable_path.as_posix()]
    elif extension == ".R":
        # Execute R script
        return ["Rscript", executable_path.as_posix()]
    elif extension == ".sh":
        # Execute shell script
        return ["sh", executable_path.as_posix()]
    else:
        logging.error(
            f"ERROR: Unsupported script extension: {extension}. Only Python/R/Shell scripts are supported."
        )
        return None
