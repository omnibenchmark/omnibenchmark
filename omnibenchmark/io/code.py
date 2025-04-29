import logging
import os
import tarfile

from pathlib import Path
from filelock import FileLock
from git import Repo
from git.exc import InvalidGitRepositoryError
import requests

from omnibenchmark.workflow.snakemake.scripts.utils import (
    generate_unique_repo_folder_name,
)


# Get the code from a GitHub repository  ( https://docs.github.com/en/rest/repos/contents?apiVersion=2022-11-28#download-a-repository-archive-tar)
def get_git_archive_from_github(
    output_dir: Path, repository_url: str, commit_hash: str
) -> Path:
    # construct url
    baseurl = (
        repository_url.replace("https://", "")
        .replace("http://", "")
        .replace("github.com/", "")
        .replace(".git", "")
    )
    url = f"https://api.github.com/repos/{baseurl}/tarball/{commit_hash}"
    headers = {
        "Accept": "application/vnd.github+json",
        "X-GitHub-Api-Version": "2022-11-28",
    }
    # get the tar file
    req = requests.get(url, headers=headers)
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
    output_dir: Path, repository_url: str, commit_hash: str
) -> Path:
    # construct url
    baseurl = repository_url.replace("https://", "").replace("http://", "")
    gitlab_remote = str(Path(baseurl).parents[-2])
    baseurl = baseurl.replace(f"{gitlab_remote}/", "")
    # url encode the baseurl
    baseurl = requests.utils.quote(baseurl, safe="")
    url = f"https://{gitlab_remote}/api/v4/projects/{str(baseurl)}/repository/archive?sha={commit_hash}"
    # get the tar file
    req = requests.get(url)
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


def get_git_archive(output_dir: Path, repository_url: str, commit_hash: str) -> Path:
    if "github.com" in repository_url:
        get_git_archive_from_github(output_dir, repository_url, commit_hash)
    elif "gitlab" in repository_url:
        get_git_archive_from_gitlab(output_dir, repository_url, commit_hash)
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
            except Exception:
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
            except Exception:
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
