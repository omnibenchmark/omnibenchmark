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
    i = 0
    while i < len(lines):
        line = lines[i]
        # Skip lines that contain [WARN]
        if "[WARN]" in line:
            i += 1
            # Skip continuation lines after [WARN] (lines that don't look like normal output)
            while (
                i < len(lines)
                and lines[i].strip()
                and not any(
                    lines[i].startswith(prefix)
                    for prefix in [
                        "Running",
                        "Found",
                        "Error:",
                        "Module",
                        "Benchmark",
                        "Building",
                        "Selecting",
                        "Job",
                        "Nothing",
                    ]
                )
            ):
                i += 1
            continue
        filtered_lines.append(line)
        i += 1
    normalized_output = "\n".join(filtered_lines)

    # Replace multiple spaces and tabs with a single space
    normalized_output = re.sub(r"[ \t]+", " ", normalized_output)

    return normalized_output
