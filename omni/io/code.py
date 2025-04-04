import hashlib
import logging
import os
import shutil
import tarfile
from pathlib import Path
from typing import List, Union
from filelock import FileLock

import requests
from git import Repo
from git.exc import InvalidGitRepositoryError
from omni.workflow.snakemake.scripts.utils import generate_unique_repo_folder_name

REPOSITORIES_DIR = Path(".snakemake") / "repos"


def construct_url_github(
    repository_url: str, commit_hash: Union[None, str] = None
) -> str:
    baseurl = (
        repository_url.replace("https://", "")
        .replace("http://", "")
        .replace("github.com/", "")
        .replace(".git", "")
    )
    if commit_hash is None:
        return f"https://api.github.com/repos/{baseurl}"
    else:
        return f"https://api.github.com/repos/{baseurl}/tarball/{commit_hash}"


def construct_url_gitlab(
    repository_url: str, commit_hash: Union[None, str] = None
) -> str:
    baseurl = repository_url.replace("https://", "").replace("http://", "")
    gitlab_remote = str(Path(baseurl).parents[-2])
    baseurl = baseurl.replace(f"{gitlab_remote}/", "")
    baseurl = requests.utils.quote(baseurl, safe="")
    if commit_hash is None:
        return f"https://{gitlab_remote}/api/v4/projects/{str(baseurl)}"
    else:
        return f"https://{gitlab_remote}/api/v4/projects/{str(baseurl)}/repository/archive?sha={commit_hash}"


def construct_url(repository_url: str, commit_hash: Union[None, str] = None) -> str:
    if "github.com" in repository_url:
        return construct_url_github(repository_url, commit_hash)
    elif "gitlab" in repository_url:
        return construct_url_gitlab(repository_url, commit_hash)
    else:
        raise ValueError(f"Unknown git provider in {repository_url}")


def check_remote_repo_existance(
    repository_url: str, commit_hash: Union[None, str] = None
) -> bool:
    url = construct_url(repository_url, commit_hash)
    req = requests.get(url)
    return req.ok


# Get the code from a GitHub repository  ( https://docs.github.com/en/rest/repos/contents?apiVersion=2022-11-28#download-a-repository-archive-tar)
def get_git_archive_from_github(
    output_dir: Path,
    repository_url: str,
    commit_hash: Union[None, str] = None,
    dry_run: bool = False,
) -> Path:
    # construct url
    url = construct_url_github(repository_url, commit_hash)
    headers = {
        "Accept": "application/vnd.github+json",
        "X-GitHub-Api-Version": "2022-11-28",
    }
    # get the tar file
    req = requests.get(url, headers=headers)
    if dry_run:
        return req.ok
    tarfilename = Path(f"{output_dir}.tar")
    with open(tarfilename, "wb") as f:
        f.write(req.content)
    # extract the tar file
    with tarfile.open(tarfilename) as f:
        members = f.getmembers()
        members_modified = []
        for member in members:
            if member.isfile():
                member.path = str(Path("/".join(Path(member.path).parts[1:])))
                members_modified.append(member)
        f.extractall(members=members_modified, path=output_dir, filter="data")
    if tarfilename.exists():
        os.remove(tarfilename)


# https://docs.gitlab.com/ee/api/repositories.html#get-file-archive
def get_git_archive_from_gitlab(
    output_dir: Path,
    repository_url: str,
    commit_hash: Union[None, str] = None,
    dry_run: bool = False,
) -> Path:
    # construct url
    url = construct_url_gitlab(repository_url, commit_hash)
    # get the tar file
    req = requests.get(url)
    if dry_run:
        return req.ok
    tarfilename = Path(f"{output_dir}.tar")
    with open(tarfilename, "wb") as f:
        f.write(req.content)
    # extract the tar file
    with tarfile.open(tarfilename) as f:
        members = f.getmembers()
        members_modified = []
        for member in members:
            if member.isfile():
                member.path = str(Path("/".join(Path(member.path).parts[1:])))
                members_modified.append(member)
        f.extractall(members=members_modified, path=output_dir, filter="data")
    if tarfilename.exists():
        os.remove(tarfilename)


def get_git_archive(
    output_dir: Path, repository_url: str, commit_hash: str, dry_run: bool = False
) -> Path:
    if "github.com" in repository_url:
        return get_git_archive_from_github(
            output_dir, repository_url, commit_hash, dry_run
        )
    elif "gitlab" in repository_url:
        return get_git_archive_from_gitlab(
            output_dir, repository_url, commit_hash, dry_run
        )
    else:
        raise ValueError(f"Unknown git provider in {repository_url}")


def clone_module(output_dir: Path, repository_url: str, commit_hash: str) -> Path:
    module_name = generate_unique_repo_folder_name(repository_url, commit_hash)
    module_dir = output_dir / module_name

    lock = module_dir.with_suffix(".lock")
    with FileLock(lock):
        if not module_dir.exists():
            try:
                logging.info(
                    f"Get git archive `{repository_url}:{commit_hash}` to `{module_dir.as_posix()}`"
                )
                get_git_archive(module_dir, repository_url, commit_hash)
                observed_commit_hash = commit_hash
            except Exception as e:
                logging.info(
                    f"Archival retirieval failed, cloning module `{repository_url}:{commit_hash}` to `{module_dir.as_posix()}`"
                )
                repo = Repo.clone_from(repository_url, module_dir.as_posix())
                repo.git.checkout(commit_hash)
                observed_commit_hash = repo.head.commit.hexsha[:7]
        else:
            try:
                repo = Repo(module_dir.as_posix())
                observed_commit_hash = repo.head.commit.hexsha[:7]
            except InvalidGitRepositoryError:
                # is archive
                observed_commit_hash = commit_hash
            except Exception as e:
                observed_commit_hash = "known"

        if observed_commit_hash != commit_hash:
            logging.error(
                f"ERROR: Failed while cloning module `{repository_url}:{commit_hash}`"
            )
            logging.error(
                f"{commit_hash} does not match {repo.head.commit.hexsha[:7]}`"
            )
            raise RuntimeError(
                f"ERROR: {commit_hash} does not match {repo.head.commit.hexsha[:7]}"
            )

    return module_dir


def dump_parameters_to_file(output_dir: str, parameters: List[str]):
    os.makedirs(output_dir, exist_ok=True)

    if parameters is not None:
        params_file = os.path.join(output_dir, "parameters.txt")
        with open(params_file, "w") as params_file:
            params_file.write(f"{parameters}")

        param_dict_file = os.path.join(output_dir, "..", "parameters_dict.txt")
        with open(param_dict_file, "a") as param_dict_file:
            param_dict_file.write(f"{os.path.basename(output_dir)} {parameters}\n")
