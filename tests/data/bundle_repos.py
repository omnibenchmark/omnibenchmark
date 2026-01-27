#!/usr/bin/env python3
"""
This script creates Git bundle files for remote repositories, so that we don't depend on the network and
we have a self-contained copy of the data.
"""

import tempfile
import subprocess
from pathlib import Path
from typing import Dict, Tuple

REPOS = {
    "https://github.com/omnibenchmark-example/data.git": "bundles/data.bundle",
    "https://github.com/omnibenchmark-example/process.git": "bundles/process.bundle",
    "https://github.com/omnibenchmark-example/method.git": "bundles/method.bundle",
    "https://github.com/omnibenchmark-example/metric.git": "bundles/metric.bundle",
    "https://github.com/btraven00/dummymodule": "bundles/dummymodule.bundle",
    "https://github.com/btraven00/ob-test-collector": "bundles/dummycollector.bundle",
    "https://github.com/btraven00/omni-module-tests": "bundles/omni-module-tests.bundle",
}


def create_single_repo_bundle(remote_url: str, bundle_basepath: str) -> Tuple[str, str]:
    """
    Create a Git bundle file for a single remote repository (HEAD only).
    The bundle filename will include the commit hash.

    Args:
        remote_url: URL of the remote repository
        bundle_basepath: Output path for the bundle (without extension - .bundle will be added)

    Returns:
        Tuple of (bundle_path, commit_hash)
    """
    bundle_basepath = Path(bundle_basepath).resolve()

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir_path = Path(temp_dir)

        # Clone the repo (shallow)
        print(f"Cloning {remote_url}...")
        subprocess.run(
            ["git", "clone", "--depth", "1", remote_url, temp_dir_path],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )

        # Get the HEAD commit hash (first 7 characters)
        commit_hash = subprocess.run(
            ["git", "rev-parse", "--short=7", "HEAD"],
            cwd=temp_dir_path,
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()

        # Create the final bundle path with commit hash
        bundle_path = (
            bundle_basepath.parent
            / f"{bundle_basepath.stem}_{commit_hash}{bundle_basepath.suffix or '.bundle'}"
        )
        bundle_path.parent.mkdir(parents=True, exist_ok=True)

        # Create the bundle
        print(f"Creating bundle: {bundle_path}")
        subprocess.run(
            ["git", "bundle", "create", bundle_path, "HEAD"],
            cwd=temp_dir_path,
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )

        print(f"Created bundle for commit: {commit_hash}\n")
        return (str(bundle_path), commit_hash)


def create_repo_bundles(repo_map: Dict[str, str]) -> Dict[str, str]:
    """
    Create Git bundle files for multiple remote repositories.

    Args:
        repo_map: Dictionary mapping remote URLs to output bundle file paths
                 (without extension - .bundle will be added)

    Returns:
        Dictionary mapping bundle file paths to their HEAD commit hashes
    """
    results = {}

    for remote_url, bundle_basepath in repo_map.items():
        try:
            bundle_path, commit_hash = create_single_repo_bundle(
                remote_url, bundle_basepath
            )
            results[bundle_path] = commit_hash
        except subprocess.CalledProcessError as e:
            print(f"Error processing {remote_url}: {e}")
            continue

    return results


if __name__ == "__main__":
    print("Creating repository bundles with commit hashes...\n")
    results = create_repo_bundles(REPOS)

    print("\nBundle creation complete. Results:")
    for bundle_path, commit_hash in results.items():
        print(f"- {Path(bundle_path).name}: {commit_hash}")
        print(
            f"  Use with: git clone {bundle_path} && cd $(basename {bundle_path} .bundle) && git checkout {commit_hash}"
        )
