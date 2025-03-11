import hashlib
import os
from pathlib import Path
from typing import List


def get_module_name_from_rule_name(rule_name: str) -> str:
    return rule_name.split("_")[1]


# Create a unique folder name based on the repository URL and commit hash
def generate_unique_repo_folder_name(repo_url: str, commit_hash: str) -> str:
    unique_string = f"{repo_url}@{commit_hash}"
    folder_name = hashlib.md5(unique_string.encode()).hexdigest()

    return folder_name


def dump_parameters_to_file(output_dir: Path, parameters: List[str]) -> None:
    os.makedirs(output_dir, exist_ok=True)

    if parameters is not None:
        params_file = os.path.join(output_dir, "parameters.txt")
        with open(params_file, "w") as params_file:
            params_file.write(f"{parameters}")

        param_dict_file = os.path.join(output_dir, "..", "parameters_dict.txt")
        with open(param_dict_file, "a") as param_dict_file:
            param_dict_file.write(f"{os.path.basename(output_dir)} {parameters}\n")
