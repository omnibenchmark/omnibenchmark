"""Error formatting for CLI output."""

from omnibenchmark.model.validation import BenchmarkParseError


def pretty_print_parse_error(error: BenchmarkParseError) -> str:
    """Format a BenchmarkParseError to present useful information to the user.

    Displays the error with file location, line marker, and context
    information.

    Args:
        error: The BenchmarkParseError to format

    Returns:
        A formatted error message with file location, line marker, and context

    Example output:
        Failed to process parameter values
            --> /path/to/file.yml:153

          153  - values: ["--method", "genie", "--threshold", 0.5]
                 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

          Stage: clustering
          Module: genieclust
          Parameter index: 0
          Values: ['--method', 'genie', '--threshold', 0.5]
          Error: '<' not supported between instances of 'float' and 'str'
    """
    message_parts = [error.message]

    # Add file location and line marker if available
    if error.yaml_file and error.line_number:
        try:
            with open(error.yaml_file, "r") as f:
                lines = f.readlines()
                if 0 < error.line_number <= len(lines):
                    # Show context: 2 lines before and 1 line after
                    context_before = 2
                    context_after = 1
                    start_line = max(1, error.line_number - context_before)
                    end_line = min(len(lines), error.line_number + context_after)

                    # Calculate padding for line numbers
                    line_num_width = len(str(end_line))

                    location_info = f"\n  --> {error.yaml_file}:{error.line_number}\n"

                    # Show lines with context
                    for line_idx in range(start_line, end_line + 1):
                        line_content = lines[line_idx - 1].rstrip()
                        location_info += (
                            f"{line_idx:>{line_num_width}} | {line_content}\n"
                        )

                        # Add marker under the offending line
                        if line_idx == error.line_number:
                            indent = len(line_content) - len(line_content.lstrip())
                            marker = " " * indent + "^" * max(
                                1, len(line_content.strip())
                            )
                            location_info += f"{' ' * line_num_width} | {marker}\n"

                    message_parts.append(location_info)
        except Exception:
            # Simple fallback if for whatever reason we cannot read the file
            message_parts.append(
                f"\n  File: {error.yaml_file}\n  Line: {error.line_number}\n"
            )

    # Add context information
    context_parts = []
    if error.stage_id:
        context_parts.append(f"  Stage: {error.stage_id}")
    if error.module_id:
        context_parts.append(f"  Module: {error.module_id}")
    if error.parameter_index is not None:
        context_parts.append(f"  Parameter index: {error.parameter_index}")
    if error.values is not None:
        context_parts.append(f"  Values: {error.values}")
    if error.original_error:
        context_parts.append(f"  Error: {error.original_error}")

    if context_parts:
        message_parts.append("\n" + "\n".join(context_parts))

    return "".join(message_parts)
