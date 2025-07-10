# Data model and functions to manage repositories

import hashlib

# TODO: model repo as its own class, for validation, housekeeping and error handling.


# Create a unique folder name based on the repository URL and commit hash
# TODO: our strategy generates unnecesary duplication of cloned repos. Could clone the repos in a cache folder instead,
# and checkout just the commit hash.
def get_repo_hash(repo: str, ref: str) -> str:
    """
    repo is the URL or local path
    ref is expected to be a commit hash or branch name
    """
    if not repo or not ref:
        raise ValueError("Both repo_url and commit_hash must be provided")
    if repo == "" or ref == "":
        raise ValueError("Both repo_url and commit_hash must be non-empty")
    if len(ref) < 7:
        raise ValueError("commit_hash must be at least 7 characters long")

    repr = f"{repo}@{ref}"
    return hashlib.sha256(repr.encode("utf-8")).hexdigest()
