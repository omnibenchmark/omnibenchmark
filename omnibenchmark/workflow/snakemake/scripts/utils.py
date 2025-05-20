import hashlib


def get_module_name_from_rule_name(rule_name: str) -> str:
    return rule_name.split("_")[1]


# Create a unique folder name based on the repository URL and commit hash
def generate_unique_repo_folder_name(repo_url: str, commit_hash: str) -> str:
    unique_string = f"{repo_url}@{commit_hash}"
    folder_name = hashlib.sha256(unique_string.encode()).hexdigest()

    return folder_name
