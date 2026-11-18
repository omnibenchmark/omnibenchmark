"""Validator for the tutorial benchmark: a scores file must not be empty.

Called with the file to check as its only argument; exit 0 means it passes."""

import sys

with open(sys.argv[1]) as f:
    lines = [line for line in f if line.strip()]

if not lines:
    sys.exit(f"{sys.argv[1]} is empty")
