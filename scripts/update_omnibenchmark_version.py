import argparse
import re
from pathlib import Path


def update_version(file_path: str, old_version: str, new_version: str):
    file = Path(file_path)
    if not file.exists():
        print(f"File not found: {file_path}")
        return

    content = file.read_text()

    # Define all replacement patterns
    patterns = [
        (rf"(omnibenchmark[=]+){re.escape(old_version)}", rf"\g<1>{new_version}"),
        (rf"(omnibenchmark[-]){re.escape(old_version)}", rf"\g<1>{new_version}"),
        (
            rf"(OmniBenchmark CLI, version )({re.escape(old_version)})",
            rf"\g<1>{new_version}",
        ),
        (
            rf"(?<![\w.-])version {re.escape(old_version)}(?![\w.-])",
            f"version {new_version}",
        ),
    ]

    updated = content
    for pattern, replacement in patterns:
        updated = re.sub(pattern, replacement, updated)

    if updated != content:
        file.write_text(updated)
        print(f"Updated version in {file_path}")
    else:
        print(f"No changes made in {file_path} (version {old_version} not found).")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Update omnibenchmark version in documentation and YAML."
    )
    parser.add_argument("old_version", help="Old version string (e.g., 0.2.1)")
    parser.add_argument("new_version", help="New version string (e.g., 0.3.0)")
    parser.add_argument(
        "files",
        nargs="+",
        help="List of files to update (e.g., howto.md omni-environment.yml)",
    )

    args = parser.parse_args()
    for f in args.files:
        update_version(f, args.old_version, args.new_version)
