"""Model-independent message formatting shared across layers."""

from pathlib import Path
from typing import Optional


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
