"""Helper functions for parsing and error handling in benchmark models."""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Any, Optional

from pydantic import ValidationError as PydanticValidationError

from .validation import BenchmarkParseError


@dataclass
class ErrorLocationIndices:
    """Indices extracted from a Pydantic error location tuple."""

    stage_idx: Optional[int] = None
    module_idx: Optional[int] = None
    param_idx: Optional[int] = None
    value_idx: Optional[int] = None


@dataclass
class ErrorContext:
    """Contextual information extracted from benchmark data."""

    stage_id: Optional[str] = None
    module_id: Optional[str] = None
    values: Optional[list] = None


def _extract_indices_from_location(loc: tuple) -> ErrorLocationIndices:
    """Extract stage, module, parameter, and value indices from Pydantic error location.

    Args:
        loc: Location tuple like ('stages', 1, 'modules', 3, 'parameters', 0, 'values', 3)

    Returns:
        ErrorLocationIndices with extracted indices
    """
    indices = ErrorLocationIndices()

    try:
        for i, part in enumerate(loc):
            if i + 1 >= len(loc):
                break
            if part == "stages":
                indices.stage_idx = loc[i + 1]
            elif part == "modules":
                indices.module_idx = loc[i + 1]
            elif part == "parameters":
                indices.param_idx = loc[i + 1]
            elif part == "values":
                indices.value_idx = loc[i + 1]
    except (IndexError, TypeError):
        pass

    return indices


def _get_context_from_data(
    data: Dict[str, Any],
    indices: ErrorLocationIndices,
) -> ErrorContext:
    """Extract stage_id, module_id, and parameter values from data.

    Args:
        data: The parsed YAML data dict
        indices: Location indices from the error

    Returns:
        ErrorContext with extracted contextual information
    """
    if indices.stage_idx is None or "stages" not in data:
        return ErrorContext()

    stage = data["stages"][indices.stage_idx]
    stage_id = stage.get("id", f"<stage index {indices.stage_idx}>")

    if indices.module_idx is None or "modules" not in stage:
        return ErrorContext(stage_id=stage_id)

    module = stage["modules"][indices.module_idx]
    module_id = module.get("id", f"<module index {indices.module_idx}>")

    if indices.param_idx is None or "parameters" not in module:
        return ErrorContext(stage_id=stage_id, module_id=module_id)

    param = module["parameters"][indices.param_idx]
    values = param.get("values")

    return ErrorContext(stage_id=stage_id, module_id=module_id, values=values)


def _find_line_number(
    line_map: Dict[str, int],
    loc: tuple,
    indices: ErrorLocationIndices,
) -> Optional[int]:
    """Find the line number in the YAML file by trying different path patterns.

    Args:
        line_map: Mapping of YAML paths to line numbers
        loc: The original error location tuple
        indices: Extracted location indices

    Returns:
        Line number if found, None otherwise
    """
    if not line_map:
        return None

    # Build path patterns from the location tuple
    patterns = []

    # For nested paths (stages/modules/parameters), build specific patterns
    if indices.stage_idx is not None and indices.param_idx is not None:
        patterns.extend(
            [
                f"stages[{indices.stage_idx}].modules[{indices.module_idx}].parameters[{indices.param_idx}].values[{indices.value_idx}]",
                f"stages[{indices.stage_idx}].modules[{indices.module_idx}].parameters[{indices.param_idx}].values",
                f"stages[{indices.stage_idx}].modules[{indices.module_idx}].parameters[{indices.param_idx}]",
            ]
        )

    # For top-level or simple paths, convert location tuple to path string
    # e.g., ('version',) -> 'version' or ('storage', 'endpoint') -> 'storage.endpoint'
    if loc:
        simple_path = ".".join(str(part) for part in loc if isinstance(part, str))
        if simple_path:
            patterns.append(simple_path)
            # For missing nested fields (e.g., storage.endpoint), also try parent path
            if "." in simple_path:
                parent_path = ".".join(simple_path.split(".")[:-1])
                patterns.append(parent_path)

    # Try all patterns
    for pattern in patterns:
        if pattern in line_map:
            return line_map[pattern]

    return None


def convert_pydantic_error_to_parse_error(
    pydantic_error: PydanticValidationError,
    data: Dict[str, Any],
    line_map: Dict[str, int],
    yaml_file_path: Optional[Path],
) -> BenchmarkParseError:
    """Convert a Pydantic validation error to a BenchmarkParseError with context.

    Args:
        pydantic_error: The Pydantic validation error
        data: The parsed YAML data dict
        line_map: Mapping of YAML paths to line numbers
        yaml_file_path: Path to the YAML file (if loaded from file)

    Returns:
        BenchmarkParseError with contextual information
    """
    if not pydantic_error.errors():
        return BenchmarkParseError(
            message=str(pydantic_error),
            yaml_file=yaml_file_path,
            original_error=pydantic_error,
        )

    error = pydantic_error.errors()[0]
    loc = error.get("loc", ())

    # Extract indices from the error location
    indices = _extract_indices_from_location(loc)

    # Get contextual information from the data
    context = _get_context_from_data(data, indices)

    # Find the line number in the YAML file
    line_num = None
    if yaml_file_path:
        line_num = _find_line_number(line_map, loc, indices)

    return BenchmarkParseError(
        message=error.get("msg", str(pydantic_error)),
        yaml_file=yaml_file_path,
        line_number=line_num,
        stage_id=context.stage_id,
        module_id=context.module_id,
        parameter_index=indices.param_idx,
        values=context.values,
        original_error=pydantic_error,
    )
