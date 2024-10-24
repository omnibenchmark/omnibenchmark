import sys

import click


## to validate the YAML
def validate_benchmark(benchmark_file: str, echo: bool = True):
    from pathlib import Path

    import yaml

    from omni.benchmark import Benchmark

    if benchmark_file.endswith(".yaml") or benchmark_file.endswith(".yml"):
        try:
            with open(benchmark_file, "r") as file:
                yaml.safe_load(file)
                benchmark = Benchmark(Path(benchmark_file))

                if echo:
                    click.echo("Benchmark YAML file integrity check passed.")

                return benchmark

        except ValueError as e:
            click.echo(
                f"Error: Failed to parse YAML as a valid OmniBenchmark: {str(e)}.",
                err=True,
            )
            sys.exit(1)  # raise click.Exit(code=1)

        except yaml.YAMLError as e:
            click.echo(f"Error: YAML file format error: {e}.", err=True)
            sys.exit(1)  # raise click.Exit(code=1)

        except FileNotFoundError:
            click.echo(
                "Error: Benchmark YAML file not found.",
                err=True,
            )
            sys.exit(1)  # raise click.Exit(code=1)

        except Exception as e:
            click.echo(
                f"Error: An unexpected error occurred: {e}",
                err=True,
            )
            sys.exit(1)  # raise click.Exit(code=1)

    else:
        click.echo(
            "Error: Invalid benchmark input. Please provide a valid YAML file path.",
            err=True,
        )
        sys.exit(1)  # raise click.Exit(code=1)
