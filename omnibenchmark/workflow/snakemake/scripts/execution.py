import logging
import os
import signal
import subprocess
import warnings
from pathlib import Path
from typing import List, Optional

import configparser
import yaml

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
    timeout: Optional[int] = None,
) -> int:
    executable = _read_entrypoint(module_dir, module_name)
    if executable is None:
        logging.error(f"ERROR: Failed to load entrypoint for {module_name}.")
        return -1  # return exit_code -1

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

    if timeout is not None:
        logging.info(f"Setting timeout to {timeout} seconds")

    exit_code = None
    process = None
    try:
        with open(stdout_file, "w") as stdout_f, open(stderr_file, "w") as stderr_f:
            # Start process with new process group to enable killing entire process tree
            process = subprocess.Popen(
                command,
                stdout=stdout_f,
                stderr=stderr_f,
                text=True,
                start_new_session=True,  # Creates new process group on POSIX
            )

            try:
                exit_code = process.wait(timeout=timeout)
                if exit_code == 0:
                    logging.info(
                        f"{executable} ran successfully with return code {exit_code}."
                    )
                else:
                    logging.error(
                        f"ERROR: Executing {executable} failed with exit code {exit_code}. "
                        f"See {stderr_file} for details."
                    )

            except subprocess.TimeoutExpired:
                timeout_msg = (
                    f"{timeout}s" if timeout is not None else "unknown timeout"
                )
                logging.error(
                    f"ERROR: Executing {executable} timed out after {timeout_msg}. "
                    f"Terminating process and all children..."
                )

                # Try graceful termination first (SIGTERM to entire process group)
                try:
                    os.killpg(os.getpgid(process.pid), signal.SIGTERM)

                    # Wait for graceful shutdown (5 seconds)
                    try:
                        process.wait(timeout=5)
                        logging.info("Process terminated gracefully.")
                    except subprocess.TimeoutExpired:
                        # Force kill if graceful termination failed
                        logging.warning(
                            "Process did not terminate gracefully, forcing kill..."
                        )
                        os.killpg(os.getpgid(process.pid), signal.SIGKILL)
                        process.wait()

                except ProcessLookupError:
                    # Process already died
                    pass

                # Return timeout error code instead of re-raising
                exit_code = (
                    124  # Standard timeout exit code (same as 'timeout' command)
                )

    finally:
        # Ensure process and all children are cleaned up
        if process is not None and process.poll() is None:
            try:
                os.killpg(os.getpgid(process.pid), signal.SIGKILL)
                process.wait()
            except (ProcessLookupError, OSError):
                pass

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


def _read_entrypoint(module_dir: Path, module_name: str) -> Optional[str]:
    """Read entrypoint from omnibenchmark.yaml or fall back to config.cfg."""
    # Try new-style omnibenchmark.yaml first
    yaml_path = module_dir / "omnibenchmark.yaml"
    if os.path.exists(yaml_path):
        try:
            with open(yaml_path, "r") as f:
                config = yaml.safe_load(f)

            if config and "entrypoints" in config:
                entrypoints = config["entrypoints"]
                if isinstance(entrypoints, dict) and "default" in entrypoints:
                    return entrypoints["default"]
                else:
                    logging.error(
                        f"ERROR: Invalid omnibenchmark.yaml format in {module_name}. "
                        "Expected 'entrypoints.default' key."
                    )
                    return None
        except yaml.YAMLError as e:
            logging.error(
                f"ERROR: Failed to parse omnibenchmark.yaml in {module_name}: {e}"
            )
            return None

    # Fall back to old-style config.cfg
    config_path = module_dir / "config.cfg"
    if os.path.exists(config_path):
        warnings.warn(
            f"Module '{module_name}' is using deprecated config.cfg. "
            "Please migrate to omnibenchmark.yaml. Support for config.cfg will be removed in a future version.",
            FutureWarning,
            stacklevel=3,
        )

        parser = configparser.ConfigParser()
        parser.read(config_path)

        if "DEFAULT" in parser and "SCRIPT" in parser["DEFAULT"]:
            return parser["DEFAULT"]["SCRIPT"]
        else:
            logging.error(f"ERROR: Invalid config.cfg format in {module_name}.")
            return None

    logging.error(
        f"ERROR: No configuration file found in {module_name}. "
        "Expected omnibenchmark.yaml or config.cfg."
    )
    return None


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
