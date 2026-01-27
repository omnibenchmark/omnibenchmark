"""Error formatting for CLI output."""

from pathlib import Path
from typing import Optional
from omnibenchmark.model.validation import BenchmarkParseError


def format_yaml_warning(
    message: str,
    yaml_file: Optional[Path] = None,
    line_number: Optional[int] = None,
    stage_id: Optional[str] = None,
    module_id: Optional[str] = None,
) -> str:
    """Format a deprecation warning with YAML context in a compact format.

    Args:
        message: The warning message
        yaml_file: Path to the YAML file
        line_number: Line number where the issue occurs
        stage_id: Optional stage ID for context
        module_id: Optional module ID for context

    Returns:
        Formatted warning string with file location and context

    Example output:
        [WARN] 05_cartesian.yaml:42 (stage 'clustering', module 'method1')
        The 'values' parameter format is deprecated.
    """
    import click

    parts = [click.style("[WARN]", fg="yellow", bold=True)]

    # Add file location
    if yaml_file:
        location = str(yaml_file)
        if line_number:
            location += f":{line_number}"
        parts.append(location)

    # Add context in compact format
    context_parts = []
    if stage_id:
        context_parts.append(f"stage '{stage_id}'")
    if module_id:
        context_parts.append(f"module '{module_id}'")

    if context_parts:
        parts.append(f"({', '.join(context_parts)})")

    # First line with location and context
    header = " ".join(parts)

    # Message on the next line
    return f"{header}\n{message}"


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
        # Clean up the error message - remove Pydantic URLs and extract key info
        error_str = str(error.original_error)
        # Remove the "For further information visit ..." line
        if "For further information visit" in error_str:
            error_str = error_str.split("For further information visit")[0].strip()

        # Try to extract just the field path and error type from Pydantic errors
        # Format: "1 validation error for Benchmark\nfield.path\n  Error message [type=...]"
        lines = error_str.strip().split("\n")
        if len(lines) >= 2 and "validation error" in lines[0]:
            # Extract field path and error message
            field_info = []
            for i in range(1, len(lines)):
                line = lines[i].strip()
                if line and not line.startswith("For further"):
                    field_info.append(line)
            if field_info:
                field_path = field_info[0]
                if len(field_info) > 1:
                    # Clean up the error description
                    detail = field_info[1]
                    # Remove the [type=...] suffix
                    if "[type=" in detail:
                        detail = detail.split("[type=")[0].strip()

                    # Make "Field required" more specific by including the field name
                    if detail.lower() == "field required":
                        # Extract just the field name from the path (e.g., "endpoint" from "storage.endpoint")
                        field_name = field_path.split(".")[-1]
                        detail = f"Missing required field '{field_name}'"

                    context_parts.append(f"  Field path: {field_path}")
                    context_parts.append(f"  Issue: {detail}")
                else:
                    context_parts.append(f"  Field: {field_path}")
        else:
            # Fallback for non-Pydantic errors
            context_parts.append(f"  Details: {error_str}")

    if context_parts:
        message_parts.append("\n" + "\n".join(context_parts))

    return "".join(message_parts)
