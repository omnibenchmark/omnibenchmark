import re

# a couple of assertion shortcuts


def assert_startswith(fd, expected):
    assert clean(fd).startswith(clean(expected))


def assert_in_output(fd, expected):
    assert clean(expected) in clean(fd)


def clean(output: str) -> str:
    output = output.strip()
    output = output.replace("    ", "")

    # Replace different newline characters with a single '\n'
    normalized_output = re.sub(r"\r\n|\r", "\n", output)

    # Replace multiple spaces and tabs with a single space
    normalized_output = re.sub(r"[ \t]+", " ", normalized_output)

    return normalized_output
