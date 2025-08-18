"""General utils functions"""

import os
import subprocess
from pathlib import Path
from typing import List, Union, Any


# Import moved to function to avoid circular import


def try_avail_envmodule(module_name: str) -> bool:
    """Check if an environment module is available on the system.

    NOTE: This function should only be used in execution contexts (e.g., BenchmarkExecution),
    not in abstract model validation. Models should remain system-agnostic.

    Args:
        module_name: Name of the environment module to check

    Returns:
        bool: True if module is available, False otherwise
    """
    env = {}
    env.update(os.environ)

    command = f"""
    . "$LMOD_PKG"/init/profile ;
    module purge ;
    module avail {module_name}"""
    result = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=True,
        text=True,
        env=env,
    )

    return "No module(s) or extension(s) found!" not in result.stderr


def as_list(input: Union[List, Any]):
    return input if isinstance(input, List) else [input]


def merge_dict_list(list_of_dicts):
    """Merge a list of dictionaries into a single dictionary."""
    merged_dict = {
        key: value for d in list_of_dicts if d is not None for key, value in d.items()
    }

    return merged_dict


def format_mc_output(output, out_dir: Path, collector_id: str):
    """Format metric collector output path.

    Args:
        output: IOFile object
        out_dir: Output directory path
        collector_id: Collector identifier
    """
    if output.path:
        o = output.path.replace("{input}", str(out_dir))
        o = o.replace("{name}", collector_id)
        return o
    else:
        return str(out_dir / output.id)
