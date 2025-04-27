import hashlib
import os
from pathlib import Path
from typing import Optional

from omni.benchmark.params import Params


def get_module_name_from_rule_name(rule_name: str) -> str:
    return rule_name.split("_")[1]


# Create a unique folder name based on the repository URL and commit hash
def generate_unique_repo_folder_name(repo_url: str, commit_hash: str) -> str:
    unique_string = f"{repo_url}@{commit_hash}"
    folder_name = hashlib.md5(unique_string.encode()).hexdigest()

    return folder_name


def create_parameters_symlink(output_dir: Path, parameters: Params) -> Optional[str]:
    os.makedirs(output_dir, exist_ok=True)

    if parameters:
        folder_identifier = parameters.generate_folder_identifier()

        real_target = output_dir
        symlink_path = output_dir.parent / folder_identifier

        if symlink_path.exists() or symlink_path.is_symlink():
            symlink_path.unlink()

        target_relative = os.path.relpath(real_target, symlink_path.parent)
        symlink_path.symlink_to(target_relative)

        return str(symlink_path)

    return None


def dump_parameters_to_file(output_dir: Path, parameters: Params) -> None:
    os.makedirs(output_dir, exist_ok=True)

    if parameters:
        params_file = output_dir / "parameters.txt"
        with open(params_file, "w") as params_file:
            params_file.write(f"{parameters}")

        param_dict_file = output_dir.parent / "parameters_dict.txt"
        basename = os.path.basename(output_dir)
        if os.path.exists(param_dict_file):
            with open(param_dict_file, "r") as file:
                existing_lines = file.readlines()

            # If the basename exists, skip writing
            if any(basename in line for line in existing_lines):
                return

        # If not found, write the new entry to parameters_dict.txt
        with open(param_dict_file, "a") as file:
            file.write(f"{basename} {parameters}\n")
