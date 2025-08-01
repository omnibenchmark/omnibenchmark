#!/usr/bin/env bash

# Exit if incorrect number of arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <old_version> <new_version>"
    echo "Example: $0 0.2.1 0.3.0"
    exit 1
fi

OLD_VERSION="$1"
NEW_VERSION="$2"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Hardcoded list of files to update
FILES=(
  "$SCRIPT_DIR/../README.md"
  "$SCRIPT_DIR/../omni-environment.yml"
)

# Add all Markdown files from docs/src
for mdfile in "$SCRIPT_DIR"/../docs/src/*.md; do
  [ -e "$mdfile" ] && FILES+=("$mdfile")
done

python "$SCRIPT_DIR/update_omnibenchmark_version.py" "$OLD_VERSION" "$NEW_VERSION" "${FILES[@]}"