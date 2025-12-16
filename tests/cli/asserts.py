import re

# a couple of assertion shortcuts


def assert_startswith(fd, expected):
    assert clean(fd).startswith(clean(expected))


def assert_in_output(fd, expected):
    assert clean(expected) in clean(fd)


def clean(output: str) -> str:
    output = output.strip()
    output = output.replace("    ", "")

    # Remove ANSI color codes (e.g., \x1b[31m, \x1b[1m, \x1b[0m)
    output = re.sub(r"\x1b\[[0-9;]*m", "", output)

    # Replace different newline characters with a single '\n'
    normalized_output = re.sub(r"\r\n|\r", "\n", output)

    # Remove warning lines (e.g., [WARN] messages and their continuation lines)
    lines = normalized_output.split("\n")
    filtered_lines = []
    skip_next = False
    for i, line in enumerate(lines):
        if skip_next:
            skip_next = False
            continue
        # Skip lines that start with [WARN] or are continuation lines after [WARN]
        if "[WARN]" in line:
            # Check if the next line is a continuation (starts with "The ")
            if i + 1 < len(lines) and lines[i + 1].strip().startswith("The "):
                skip_next = True
            continue
        filtered_lines.append(line)
    normalized_output = "\n".join(filtered_lines)

    # Replace multiple spaces and tabs with a single space
    normalized_output = re.sub(r"[ \t]+", " ", normalized_output)

    return normalized_output
