# validation_engine.py
#
# A module to validate artifacts produced by omnibenchmark modules.
# Internally uses the `returns` library for robust error handling.
#
# Dependencies:
# - pyyaml (already in project deps)
# - returns (only requires typing-extensions)

import csv
import hashlib
import json
import re
import yaml
from pathlib import Path
from typing import Any, Callable, TypeAlias, TypedDict
from dataclasses import dataclass, field

from returns.result import Result, Success, Failure
from omnibenchmark.core.interfaces import ValidatableModule, ValidatableBenchmark


# First, a few custom exceptions to allow us defining a hierarchy of validation errors


class ValidationError(Exception):
    """Base exception for all validation errors."""

    pass


class ArtifactNotFoundError(ValidationError):
    """Raised when an artifact file or directory does not exist."""

    pass


class RuleTypeError(ValidationError):
    """Raised when a validation rule is of an unknown type."""

    pass


class FileSystemError(ValidationError):
    """Raised for file system mismatches (e.g., expected file, got dir)."""

    pass


class ContentError(ValidationError):
    """Raised for errors related to file content (e.g., regex mismatch, checksum fail)."""

    pass


class ParsingError(ValidationError):
    """Raised for errors during file parsing (e.g., JSON, YAML, CSV)."""

    pass


# Validation Rules as type definitions


class IsFileRule(TypedDict):
    is_file: bool


class IsDirRule(TypedDict):
    is_dir: bool


class MatchesRule(TypedDict):
    matches: str


class ChecksumRule(TypedDict):
    checksum: str


class HasColumnRule(TypedDict):
    has_column: str


class NotNARule(TypedDict):
    not_na: str


class HasKeyRule(TypedDict):
    has_key: str


ValidationRule: TypeAlias = (
    IsFileRule
    | IsDirRule
    | MatchesRule
    | ChecksumRule
    | HasColumnRule
    | NotNARule
    | HasKeyRule
)
ValidationFileContent: TypeAlias = dict[str, list[ValidationRule]]


# We will use data classes for Results


@dataclass(frozen=True, slots=True)
class ValidationResult:
    """Represents the outcome of a single validation check."""

    artifact_name: str
    rule: ValidationRule
    success: bool
    message: str = ""


@dataclass(frozen=True, slots=True)
class ModuleValidationResult:
    """Aggregates all validation results for a single module."""

    module_name: str
    results: list[ValidationResult] = field(default_factory=list)
    error: str | None = None

    @property
    def all_successful(self) -> bool:
        return self.error is None and all(r.success for r in self.results)


# Individual validation handlers
# Each handler returns a Result[SuccessMessage, ValidationError]

HandlerResult: TypeAlias = Result[str, ValidationError]


def _handle_is_file(artifact_path: Path, **_: Any) -> HandlerResult:
    if not artifact_path.exists():
        return Failure(
            ArtifactNotFoundError(f"Artifact does not exist at {artifact_path}")
        )
    return (
        Success("Artifact is a file.")
        if artifact_path.is_file()
        else Failure(FileSystemError(f"Artifact is not a file: {artifact_path}"))
    )


def _handle_is_dir(artifact_path: Path, **_: Any) -> HandlerResult:
    if not artifact_path.exists():
        return Failure(
            ArtifactNotFoundError(f"Artifact does not exist at {artifact_path}")
        )
    return (
        Success("Artifact is a directory.")
        if artifact_path.is_dir()
        else Failure(FileSystemError(f"Artifact is not a directory: {artifact_path}"))
    )


def _handle_matches(artifact_path: Path, rule_value: str, **_: Any) -> HandlerResult:
    try:
        content = artifact_path.read_text()
        if re.search(rule_value, content):
            return Success("Content matches pattern.")
        return Failure(ContentError(f"Content does not match pattern: '{rule_value}'"))
    except Exception as e:
        return Failure(ContentError(f"Failed to read or match content: {e}"))


def _handle_checksum(artifact_path: Path, rule_value: str, **_: Any) -> HandlerResult:
    try:
        algo, expected_hash = rule_value.split(":", 1)
        hasher = hashlib.new(algo)
        with open(artifact_path, "rb") as f:
            while chunk := f.read(8192):
                hasher.update(chunk)
        actual_hash = hasher.hexdigest()
        if actual_hash == expected_hash:
            return Success(f"Checksum ({algo}) matches.")
        return Failure(
            ContentError(
                f"Checksum mismatch. Expected: {expected_hash}, Got: {actual_hash}"
            )
        )
    except Exception as e:
        return Failure(ContentError(f"Failed to compute checksum: {e}"))


def _handle_has_column(artifact_path: Path, rule_value: str, **_: Any) -> HandlerResult:
    try:
        with open(artifact_path, "r", newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            fieldnames = reader.fieldnames
            if not fieldnames:
                return Failure(
                    ContentError(f"CSV file has no headers: {artifact_path}")
                )
            if rule_value in fieldnames:
                return Success(f"Column '{rule_value}' found.")
        return Failure(
            ContentError(f"Column '{rule_value}' not found in {artifact_path}")
        )
    except Exception as e:
        return Failure(ParsingError(f"Failed to parse CSV: {e}"))


def _handle_not_na(artifact_path: Path, rule_value: str, **_: Any) -> HandlerResult:
    try:
        with open(artifact_path, "r", newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            fieldnames = reader.fieldnames
            if not fieldnames or rule_value not in fieldnames:
                return Failure(
                    ContentError(
                        f"Column '{rule_value}' not found, cannot check for N/A values."
                    )
                )

            # Check each row for empty/NA values in the specified column
            for row_num, row in enumerate(
                reader, start=2
            ):  # start=2 because row 1 is headers
                cell_value = row.get(rule_value, "")
                # Consider empty string, None, or common NA representations as N/A
                if not cell_value or cell_value.strip().upper() in (
                    "NA",
                    "N/A",
                    "NULL",
                    "NONE",
                ):
                    return Failure(
                        ContentError(
                            f"Column '{rule_value}' contains N/A values at row {row_num}."
                        )
                    )

        return Success(f"Column '{rule_value}' has no N/A values.")
    except Exception as e:
        return Failure(ParsingError(f"Failed to parse CSV: {e}"))


def _handle_has_key(artifact_path: Path, rule_value: str, **_: Any) -> HandlerResult:
    try:
        with artifact_path.open("r") as f:
            data = json.load(f)
        if not isinstance(data, dict):
            return Failure(ParsingError("JSON content is not an object."))
        if rule_value in data:
            return Success(f"Key '{rule_value}' found.")
        return Failure(ContentError(f"Key '{rule_value}' not found in JSON object."))
    except Exception as e:
        return Failure(ParsingError(f"Failed to parse JSON: {e}"))


# Validator Registry and function call dispatcher
# In the future, we can extend the registry at runtime with more validators, even
# validators that are not known to the package. That will nee a higher level of trust in the remote code,
# since contributed modules may contain malicious code.

VALIDATOR_REGISTRY: dict[str, Callable[..., HandlerResult]] = {
    "is_file": _handle_is_file,
    "is_dir": _handle_is_dir,
    "matches": _handle_matches,
    "checksum": _handle_checksum,
    "has_column": _handle_has_column,
    "not_na": _handle_not_na,
    "has_key": _handle_has_key,
}


def apply_rule(
    artifact_name: str, artifact_path: Path, rule: ValidationRule, module_path: Path
) -> ValidationResult:
    """Applies a single validation rule and unwraps the Result into a ValidationResult."""
    rule_type = next(iter(rule))
    # TODO: TypedDict unions don't support this pattern well at type-check time
    # Consider runtime rule validation or restructure if type safety becomes critical
    rule_value = rule[rule_type]  # type: ignore

    handler = VALIDATOR_REGISTRY.get(rule_type)
    if not handler:
        result: HandlerResult = Failure(
            RuleTypeError(f"Unknown rule type: '{rule_type}'")
        )
    else:
        result = handler(
            artifact_path=artifact_path, rule_value=rule_value, module_path=module_path
        )

    # Unwrap the result into our data class
    # Note: @safe decorator may double-wrap Success results, so we need to flatten
    if isinstance(result, Success):
        inner = result.unwrap()
        # If @safe wrapped a Success(msg), unwrap it
        if isinstance(inner, Success):
            return ValidationResult(artifact_name, rule, True, inner.unwrap())
        # If @safe wrapped a Failure, unwrap it
        elif isinstance(inner, Failure):
            return ValidationResult(artifact_name, rule, False, str(inner.failure()))
        # Normal Success case (no @safe or non-Result return)
        else:
            return ValidationResult(artifact_name, rule, True, str(inner))
    elif isinstance(result, Failure):
        return ValidationResult(artifact_name, rule, False, str(result.failure()))
    else:
        return ValidationResult(artifact_name, rule, False, "Unknown error")


def load_validation_rules(
    module_path: Path,
) -> Result[ValidationFileContent, ValidationError]:
    """Loads and parses the validation.yaml file for a module, returning a Result.

    Note: This function manually constructs Result types and handles exceptions
    internally, so @safe decorator is not needed and would cause double-wrapping.
    """
    validation_file = module_path / "validation.yaml"
    if not validation_file.is_file():
        return Failure(ArtifactNotFoundError("validation.yaml not found."))

    try:
        with validation_file.open("r") as f:
            content: ValidationFileContent = yaml.safe_load(f)
            return Success(content)
    except Exception as e:
        return Failure(ParsingError(f"Failed to parse validation.yaml: {e}"))


def run_validation_for_module(module: ValidatableModule) -> ModuleValidationResult:
    """Runs all validation checks for a single module."""
    rules_result = load_validation_rules(module.path)

    if isinstance(rules_result, Failure):
        return ModuleValidationResult(
            module_name=module.name, error=str(rules_result.failure())
        )
    elif isinstance(rules_result, Success):
        rules_content = rules_result.unwrap()
        if not rules_content:
            return ModuleValidationResult(module_name=module.name)

        all_results: list[ValidationResult] = []
        for artifact_name, rules in rules_content.items():
            artifact_path = module.artifact_paths.get(artifact_name)
            if not artifact_path:
                all_results.append(
                    ValidationResult(
                        artifact_name,
                        {"error": "Artifact not found"},  # type: ignore[typeddict-item]
                        False,
                        f"Artifact '{artifact_name}' not found in module outputs.",
                    )
                )
                continue
            for rule in rules:
                all_results.append(
                    apply_rule(artifact_name, artifact_path, rule, module.path)
                )

        return ModuleValidationResult(module.name, results=all_results)
    else:
        return ModuleValidationResult(module_name=module.name, error="Unknown error")


def run_engine(benchmark: ValidatableBenchmark) -> list[ModuleValidationResult]:
    """Main entry point to run validation across all modules in a benchmark."""
    return [run_validation_for_module(module) for module in benchmark]
