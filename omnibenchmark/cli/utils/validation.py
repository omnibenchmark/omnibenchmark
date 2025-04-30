import sys

import click
import yaml

from pathlib import Path

from omnibenchmark.benchmark import Benchmark
from omnibenchmark.cli.utils.logging import logger


## to validate the YAML
def validate_benchmark(benchmark_file: str, echo: bool = True) -> Benchmark | None:
    if benchmark_file.endswith(".yaml") or benchmark_file.endswith(".yml"):
        try:
            with open(benchmark_file, "r") as file:
                yaml.safe_load(file)
                benchmark = Benchmark(Path(benchmark_file))

                if echo:
                    logger.info("Benchmark YAML file integrity check passed.")

                return benchmark

        except ValueError as e:
            logger.error(
                f"Error: Failed to parse YAML as a valid OmniBenchmark: {str(e)}.",
            )
            sys.exit(1)

        except yaml.YAMLError as e:
            logger.info(f"Error: YAML file format error: {e}.")
            click.Abort()

        except FileNotFoundError:
            logger.error(
                "Error: Benchmark YAML file not found.",
            )
            sys.exit(1)

        except Exception as e:
            logger.error(
                f"Error: An unexpected error occurred: {e}",
            )
            sys.exit(1)

    else:
        logger.error(
            "Error: Invalid benchmark input. Please provide a valid YAML file path.",
        )
        sys.exit(1)
