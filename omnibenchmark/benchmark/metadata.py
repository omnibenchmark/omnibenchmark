"""Core validation system with strategy pattern for handling errors vs warnings.

This module provides a unified validation system that can operate in strict mode
(fail fast on errors) or warning mode (convert errors to warnings and continue).
"""

from abc import ABC, abstractmethod
from enum import Enum
from typing import Dict, List, Set, Optional, Any, Tuple
from pathlib import Path
import yaml
from spdx_license_list import LICENSES


class ValidationSeverity(Enum):
    """Severity levels for validation issues."""

    ERROR = "error"
    WARNING = "warning"


class ValidationIssue:
    """Unified validation issue that can be either an error or warning."""

    def __init__(
        self,
        issue_type: str,
        severity: ValidationSeverity,
        path: str,
        module_id: str,
        message: str,
        **context,
    ):
        self.issue_type = issue_type
        self.severity = severity
        self.path = path
        self.module_id = module_id
        self.message = message
        self.context = context

    @property
    def is_error(self) -> bool:
        """Check if this issue is an error."""
        return self.severity == ValidationSeverity.ERROR

    @property
    def is_warning(self) -> bool:
        """Check if this issue is a warning."""
        return self.severity == ValidationSeverity.WARNING

    @property
    def msg(self) -> str:
        """User-facing message (backward compatibility)."""
        return self.message

    def __eq__(self, other) -> bool:
        """Issues are equal if they have same type, path, module, and context."""
        if not isinstance(other, ValidationIssue):
            return False
        return (
            self.issue_type == other.issue_type
            and self.path == other.path
            and self.module_id == other.module_id
            and self.context == other.context
        )

    def __hash__(self) -> int:
        """Hash based on type, path, module, and context for deduplication."""
        context_items = tuple(sorted(self.context.items()))
        return hash((self.issue_type, self.path, self.module_id, context_items))

    def __repr__(self) -> str:
        return f"ValidationIssue({self.issue_type}, {self.severity.value}, {self.path}, {self.module_id})"


class ValidationStrategy(ABC):
    """Strategy for handling validation issues."""

    @abstractmethod
    def handle_issue(
        self, issue_type: str, path: str, module_id: str, message: str, **context
    ) -> ValidationIssue:
        """Handle a validation issue by creating appropriate ValidationIssue."""
        pass

    @abstractmethod
    def should_fail_fast(self) -> bool:
        """Whether to fail fast on first error."""
        pass


class StrictValidationStrategy(ValidationStrategy):
    """Strategy that maintains natural error/warning semantics but fails fast on errors."""

    def handle_issue(
        self, issue_type: str, path: str, module_id: str, message: str, **context
    ) -> ValidationIssue:
        # Determine natural severity based on issue type
        severity = self._get_natural_severity(issue_type)
        return ValidationIssue(
            issue_type=issue_type,
            severity=severity,
            path=path,
            module_id=module_id,
            message=message,
            **context,
        )

    def should_fail_fast(self) -> bool:
        return True

    def _get_natural_severity(self, issue_type: str) -> ValidationSeverity:
        """Determine the natural severity for an issue type."""
        # These are naturally warnings regardless of strategy
        warning_types = {
            "citation_author_missing_given_name",
            "omnibenchmark_yaml_invalid_yaml",
            "omnibenchmark_yaml_not_object",
            "license_file_no_license_language",
            "license_in_file_but_not_citation",
            "license_in_citation_but_not_file",
            "license_content_mismatch",
            "no_license_file",
        }

        if issue_type in warning_types:
            return ValidationSeverity.WARNING
        else:
            return ValidationSeverity.ERROR


class WarnValidationStrategy(ValidationStrategy):
    """Strategy that converts errors to warnings and continues processing."""

    def handle_issue(
        self, issue_type: str, path: str, module_id: str, message: str, **context
    ) -> ValidationIssue:
        # Certain issues are ALWAYS errors, even in warn mode
        always_error_types = {
            "omnibenchmark_yaml_missing",
        }

        if issue_type in always_error_types:
            severity = ValidationSeverity.ERROR
        else:
            # Convert other issues to warnings in warn mode
            severity = ValidationSeverity.WARNING

        return ValidationIssue(
            issue_type=issue_type,
            severity=severity,
            path=path,
            module_id=module_id,
            message=message,
            **context,
        )

    def should_fail_fast(self) -> bool:
        return False


class ValidationException(Exception):
    """Exception raised when validation fails in strict mode."""

    def __init__(self, message: str, issues: List[ValidationIssue]):
        super().__init__(message)
        self.issues = issues


class ValidationResult:
    """Container for validation results with deduplication and aggregation."""

    def __init__(self):
        self._issues: Set[ValidationIssue] = set()

    def add_issue(self, issue: ValidationIssue):
        """Add a validation issue with automatic deduplication."""
        self._issues.add(issue)

    def add_issues(self, issues: List[ValidationIssue]):
        """Add multiple validation issues."""
        for issue in issues:
            self.add_issue(issue)

    @property
    def issues(self) -> List[ValidationIssue]:
        """Get all issues sorted by module, then by type."""
        return sorted(self._issues, key=lambda x: (x.module_id, x.issue_type, x.path))

    @property
    def errors(self) -> List[ValidationIssue]:
        """Get only error issues."""
        return [issue for issue in self.issues if issue.is_error]

    @property
    def warnings(self) -> List[ValidationIssue]:
        """Get only warning issues."""
        return [issue for issue in self.issues if issue.is_warning]

    def is_valid(self) -> bool:
        """Return True if no errors."""
        return len(self.errors) == 0

    def has_warnings(self) -> bool:
        """Return True if warnings exist."""
        return len(self.warnings) > 0

    def get_issues_by_module(self) -> Dict[str, List[ValidationIssue]]:
        """Group issues by module ID."""
        result = {}
        for issue in self.issues:
            module_id = issue.module_id
            if module_id not in result:
                result[module_id] = []
            result[module_id].append(issue)
        return result

    def get_error_messages(self) -> List[str]:
        """Get user-facing error messages (backward compatibility)."""
        return [error.message for error in self.errors]

    def get_warning_messages(self) -> List[str]:
        """Get user-facing warning messages (backward compatibility)."""
        return [warning.message for warning in self.warnings]

    def add_error(
        self,
        message: str,
        issue_type: str = "generic_error",
        path: str = ".",
        module_id: str = "unknown",
    ):
        """Convenience method to add an error issue."""
        error_issue = ValidationIssue(
            issue_type=issue_type,
            severity=ValidationSeverity.ERROR,
            path=path,
            module_id=module_id,
            message=message,
        )
        self.add_issue(error_issue)

    def add_warning(
        self,
        message: str,
        issue_type: str = "generic_warning",
        path: str = ".",
        module_id: str = "unknown",
    ):
        """Convenience method to add a warning issue."""
        warning_issue = ValidationIssue(
            issue_type=issue_type,
            severity=ValidationSeverity.WARNING,
            path=path,
            module_id=module_id,
            message=message,
        )
        self.add_issue(warning_issue)


class ValidationContext:
    """Context for validation operations with strategy pattern."""

    def __init__(
        self, strategy: ValidationStrategy, module_id: str, base_path: str = ""
    ):
        self.strategy = strategy
        self.module_id = module_id
        self.base_path = base_path
        self.result = ValidationResult()

    @property
    def warn_mode(self) -> bool:
        """Check if context is in warn mode."""
        return isinstance(self.strategy, WarnValidationStrategy)

    def add_issue(self, issue_type: str, path: str, message: str, **context):
        """Add a validation issue using the current strategy."""
        full_path = str(Path(self.base_path) / path) if self.base_path else path
        issue = self.strategy.handle_issue(
            issue_type=issue_type,
            path=full_path,
            module_id=self.module_id,
            message=message,
            **context,
        )
        self.result.add_issue(issue)

        # Fail fast if strategy requires it and we have an error
        if self.strategy.should_fail_fast() and issue.is_error:
            raise ValidationException(f"Validation failed: {message}", [issue])

    def merge_result(self, other_result: ValidationResult):
        """Merge another validation result into this context."""
        self.result.add_issues(other_result.issues)


# Validation Functions
def validate_citation_cff_content(
    content: str, ctx: ValidationContext, file_path: str = "CITATION.cff"
) -> Tuple[ValidationResult, Optional[Dict[str, Any]]]:
    """Validate CITATION.cff file content."""
    if not content or not content.strip():
        ctx.add_issue(
            issue_type="citation_missing",
            path=file_path,
            message="CITATION.cff file is empty",
        )
        return ctx.result, None

    try:
        data = yaml.safe_load(content)
    except yaml.YAMLError as e:
        ctx.add_issue(
            issue_type="citation_invalid_yaml",
            path=file_path,
            message=f"Invalid YAML syntax in CITATION.cff: {str(e)}",
            yaml_error=str(e),
        )
        return ctx.result, None

    # Validate the parsed data
    validation_result, validated_data = validate_citation_cff_data(data, ctx, file_path)
    ctx.merge_result(validation_result)

    return ctx.result, validated_data


def validate_citation_cff_data(
    data: Any, ctx: ValidationContext, file_path: str = "CITATION.cff"
) -> Tuple[ValidationResult, Optional[Dict[str, Any]]]:
    """Validate parsed CITATION.cff data structure."""
    result = ValidationResult()

    if not isinstance(data, dict):
        ctx.add_issue(
            issue_type="citation_not_object",
            path=file_path,
            message="CITATION.cff content must be a YAML object/dictionary",
        )
        return result, None

    # Check required fields
    required_fields = ["cff-version", "message", "title", "authors"]
    for field in required_fields:
        if field not in data:
            ctx.add_issue(
                issue_type="citation_missing_required_field",
                path=file_path,
                message=f"Required field '{field}' missing from CITATION.cff",
                field=field,
            )

    # Validate CFF version
    if "cff-version" in data:
        version = data["cff-version"]
        if not isinstance(version, (str, float)) or str(version) not in [
            "1.0.3",
            "1.1.0",
            "1.2.0",
        ]:
            ctx.add_issue(
                issue_type="citation_unsupported_version",
                path=file_path,
                message=f"Unsupported CFF version: {version}. Supported versions: 1.0.3, 1.1.0, 1.2.0",
                version=str(version),
            )

    # Validate authors
    if "authors" in data:
        authors_result = validate_citation_authors(data["authors"], ctx, file_path)
        result.add_issues(authors_result.issues)

    # Validate license
    if "license" in data:
        license_result = validate_citation_license(data["license"], ctx, file_path)
        result.add_issues(license_result.issues)

    return result, data


def validate_citation_authors(
    authors: Any, ctx: ValidationContext, file_path: str = "CITATION.cff"
) -> ValidationResult:
    """Validate authors field in CITATION.cff."""
    result = ValidationResult()

    if authors is None:
        ctx.add_issue(
            issue_type="citation_missing_authors",
            path=file_path,
            message="Authors field is missing or null in CITATION.cff",
        )
        return result

    if not isinstance(authors, list):
        ctx.add_issue(
            issue_type="citation_authors_not_list",
            path=file_path,
            message="Authors field must be a list in CITATION.cff",
        )
        return result

    if len(authors) == 0:
        ctx.add_issue(
            issue_type="citation_missing_authors",
            path=file_path,
            message="Authors list is empty in CITATION.cff",
        )
        return result

    # Validate each author
    for i, author in enumerate(authors):
        author_result = validate_single_author(author, ctx, file_path, i)
        result.add_issues(author_result.issues)

    return result


def validate_single_author(
    author: Any,
    ctx: ValidationContext,
    file_path: str = "CITATION.cff",
    author_index: int = 0,
) -> ValidationResult:
    """Validate a single author entry."""
    result = ValidationResult()

    if not isinstance(author, dict):
        ctx.add_issue(
            issue_type="citation_author_not_object",
            path=file_path,
            message=f"Author at index {author_index} must be an object/dictionary",
            author_index=author_index,
        )
        return result

    # Check for family-names (required)
    if "family-names" not in author:
        ctx.add_issue(
            issue_type="citation_author_missing_family_name",
            path=file_path,
            message=f"Author at index {author_index} is missing required 'family-names' field",
            author_index=author_index,
        )

    # Check for given-names (warning only)
    if "given-names" not in author:
        # This should be a warning regardless of strategy, so we create it directly
        warning_issue = ValidationIssue(
            issue_type="citation_author_missing_given_name",
            severity=ValidationSeverity.WARNING,
            path=file_path,
            module_id=ctx.module_id,
            message=f"Author at index {author_index} is missing 'given-names' field (recommended)",
            author_index=author_index,
        )
        result.add_issue(warning_issue)

    return result


def validate_citation_license(
    license_value: Any, ctx: ValidationContext, file_path: str = "CITATION.cff"
) -> ValidationResult:
    """Validate license field in CITATION.cff."""
    result = ValidationResult()

    if license_value is None:
        ctx.add_issue(
            issue_type="citation_missing_license",
            path=file_path,
            message="License field is missing or null in CITATION.cff",
        )
        return result

    if not isinstance(license_value, str) or not license_value.strip():
        ctx.add_issue(
            issue_type="citation_missing_license",
            path=file_path,
            message="License field must be a non-empty string in CITATION.cff",
        )
        return result

    # Validate SPDX license identifier
    license_id = license_value.strip()
    if license_id not in LICENSES:
        ctx.add_issue(
            issue_type="citation_invalid_license",
            path=file_path,
            message=f"Invalid SPDX license identifier: {license_id}",
            license_id=license_id,
        )

    return result


def validate_omnibenchmark_yaml_content(
    content: str, ctx: ValidationContext, file_path: str = "omnibenchmark.yaml"
) -> Tuple[ValidationResult, Optional[Dict[str, Any]]]:
    """Validate omnibenchmark.yaml file content."""
    result = ValidationResult()

    if not content or not content.strip():
        # This is typically a warning, not an error
        warning_issue = ValidationIssue(
            issue_type="omnibenchmark_yaml_missing",
            severity=ValidationSeverity.WARNING,
            path=file_path,
            module_id=ctx.module_id,
            message="omnibenchmark.yaml file is empty or missing",
        )
        result.add_issue(warning_issue)
        ctx.merge_result(result)
        return result, None

    try:
        data = yaml.safe_load(content)
    except yaml.YAMLError as e:
        # YAML parsing errors are typically warnings for omnibenchmark.yaml
        warning_issue = ValidationIssue(
            issue_type="omnibenchmark_yaml_invalid_yaml",
            severity=ValidationSeverity.WARNING,
            path=file_path,
            module_id=ctx.module_id,
            message=f"Invalid YAML syntax in omnibenchmark.yaml: {str(e)}",
            yaml_error=str(e),
        )
        result.add_issue(warning_issue)
        ctx.merge_result(result)
        return result, None

    # Validate the parsed data
    validation_result, validated_data = validate_omnibenchmark_yaml_data(
        data, ctx, file_path
    )
    result.add_issues(validation_result.issues)

    ctx.merge_result(result)
    return result, validated_data


def validate_omnibenchmark_yaml_data(
    data: Any, ctx: ValidationContext, file_path: str = "omnibenchmark.yaml"
) -> Tuple[ValidationResult, Optional[Dict[str, Any]]]:
    """Validate parsed omnibenchmark.yaml data structure."""
    result = ValidationResult()

    if not isinstance(data, dict):
        # This is typically a warning
        warning_issue = ValidationIssue(
            issue_type="omnibenchmark_yaml_not_object",
            severity=ValidationSeverity.WARNING,
            path=file_path,
            module_id=ctx.module_id,
            message="omnibenchmark.yaml content must be a YAML object/dictionary",
        )
        result.add_issue(warning_issue)
        return result, None

    # Check for required fields - entrypoints with default
    if "entrypoints" not in data:
        warning_issue = ValidationIssue(
            issue_type="omnibenchmark_yaml_missing_entrypoints",
            severity=ValidationSeverity.WARNING,
            path=file_path,
            module_id=ctx.module_id,
            message="omnibenchmark.yaml is missing required 'entrypoints' field",
        )
        result.add_issue(warning_issue)
    elif not isinstance(data["entrypoints"], dict):
        warning_issue = ValidationIssue(
            issue_type="omnibenchmark_yaml_invalid_entrypoints",
            severity=ValidationSeverity.WARNING,
            path=file_path,
            module_id=ctx.module_id,
            message="omnibenchmark.yaml 'entrypoints' field must be a dictionary",
        )
        result.add_issue(warning_issue)
    elif "default" not in data["entrypoints"]:
        warning_issue = ValidationIssue(
            issue_type="omnibenchmark_yaml_missing_default_entrypoint",
            severity=ValidationSeverity.WARNING,
            path=file_path,
            module_id=ctx.module_id,
            message="omnibenchmark.yaml 'entrypoints' is missing required 'default' key",
        )
        result.add_issue(warning_issue)

    # Validate version format (should be a string)
    if "version" in data and not isinstance(data["version"], str):
        warning_issue = ValidationIssue(
            issue_type="omnibenchmark_yaml_invalid_version_format",
            severity=ValidationSeverity.WARNING,
            path=file_path,
            module_id=ctx.module_id,
            message=f"omnibenchmark.yaml 'version' field should be a string, got {type(data['version']).__name__}",
        )
        result.add_issue(warning_issue)

    return result, data


def validate_license_file_content(
    content: str, ctx: ValidationContext, file_path: str = "LICENSE"
) -> ValidationResult:
    """Validate LICENSE file content."""
    result = ValidationResult()

    if not content or not content.strip():
        ctx.add_issue(
            issue_type="license_file_empty",
            path=file_path,
            message="LICENSE file is empty",
        )
        return result

    # Check if content contains license-related keywords
    license_keywords = ["license", "copyright", "permission", "warranty", "liability"]
    content_lower = content.lower()

    if not any(keyword in content_lower for keyword in license_keywords):
        # This is typically a warning
        warning_issue = ValidationIssue(
            issue_type="license_file_no_license_language",
            severity=ValidationSeverity.WARNING,
            path=file_path,
            module_id=ctx.module_id,
            message="LICENSE file does not contain typical license language",
        )
        result.add_issue(warning_issue)

    ctx.merge_result(result)
    return result


def validate_license_consistency(
    citation_license: Optional[str],
    license_content: Optional[str],
    ctx: ValidationContext,
) -> ValidationResult:
    """Validate consistency between CITATION.cff license and LICENSE file."""
    result = ValidationResult()

    if not citation_license and not license_content:
        ctx.add_issue(
            issue_type="no_license_information",
            path=".",
            message="No license information found in either CITATION.cff or LICENSE file",
        )
        return result

    # If no citation_license but LICENSE file exists, warn about missing license in citation
    if not citation_license and license_content:
        warning_issue = ValidationIssue(
            issue_type="license_in_file_but_not_citation",
            severity=ValidationSeverity.WARNING,
            path="CITATION.cff",
            module_id=ctx.module_id,
            message="LICENSE file exists but no license specified in CITATION.cff",
        )
        result.add_issue(warning_issue)
        return result

    # If citation_license exists but no LICENSE file, that's OK - don't warn
    if citation_license and not license_content:
        return result

    # Both exist - check for consistency
    if citation_license and license_content:
        detected_license = detect_license_from_content(license_content)
        if detected_license and detected_license.upper() != citation_license.upper():
            # This is typically a warning
            warning_issue = ValidationIssue(
                issue_type="license_content_mismatch",
                severity=ValidationSeverity.WARNING,
                path=".",
                module_id=ctx.module_id,
                message=f"License mismatch: CITATION.cff specifies '{citation_license}' but LICENSE file appears to be '{detected_license}'",
                citation_license=citation_license,
                detected_license=detected_license,
            )
            result.add_issue(warning_issue)

    return result


def detect_license_from_content(content: str) -> Optional[str]:
    """Detect license type from LICENSE file content."""
    if not content:
        return None

    content_lower = content.lower()

    # Simple keyword-based detection
    if "mit license" in content_lower or "mit" in content_lower:
        return "MIT"
    elif "apache license" in content_lower or "apache" in content_lower:
        return "Apache-2.0"
    elif "gnu general public license" in content_lower or "gpl" in content_lower:
        if "version 3" in content_lower or "v3" in content_lower:
            return "GPL-3.0"
        elif "version 2" in content_lower or "v2" in content_lower:
            return "GPL-2.0"
        else:
            return "GPL"
    elif "bsd license" in content_lower or "bsd" in content_lower:
        return "BSD"

    return None


def validate_file_structure(
    files_present: Dict[str, bool],
    ctx: ValidationContext,
    citation_has_license: bool = False,
) -> ValidationResult:
    """Validate expected file structure."""
    result = ValidationResult()

    # Check for CITATION.cff
    if not files_present.get("CITATION.cff", False):
        ctx.add_issue(
            issue_type="citation_missing",
            path="CITATION.cff",
            message="CITATION.cff file not found",
        )

    # Check for LICENSE file (various names) - only warn if also missing from CITATION.cff
    license_files = [
        "LICENSE",
        "LICENSE.txt",
        "LICENSE.md",
        "LICENCE",
        "LICENCE.txt",
        "LICENCE.md",
    ]
    license_found = any(files_present.get(name, False) for name in license_files)

    # Only warn about missing LICENSE file if there's also no license in CITATION.cff
    if not license_found and not citation_has_license:
        warning_issue = ValidationIssue(
            issue_type="no_license_file",
            severity=ValidationSeverity.WARNING,
            path="LICENSE",
            module_id=ctx.module_id,
            message="No LICENSE file found",
        )
        result.add_issue(warning_issue)

    # Check for omnibenchmark.yaml or legacy config.cfg
    if not files_present.get("omnibenchmark.yaml", False):
        # Check for legacy config.cfg
        if files_present.get("config.cfg", False):
            warning_issue = ValidationIssue(
                issue_type="omnibenchmark_legacy_config",
                severity=ValidationSeverity.WARNING,
                path="config.cfg",
                module_id=ctx.module_id,
                message="Using legacy config.cfg file. Please migrate to omnibenchmark.yaml",
            )
            result.add_issue(warning_issue)
        else:
            # Missing omnibenchmark.yaml is an error - modules must have this file
            ctx.add_issue(
                issue_type="omnibenchmark_yaml_missing",
                path="omnibenchmark.yaml",
                message="omnibenchmark.yaml file not found",
            )

    return result


# Factory function for creating validation contexts
def create_validation_context(
    module_id: str, warn_mode: bool = False, base_path: str = ""
) -> ValidationContext:
    """Create a validation context with appropriate strategy."""
    if warn_mode:
        strategy = WarnValidationStrategy()
    else:
        strategy = StrictValidationStrategy()

    return ValidationContext(strategy, module_id, base_path)


# Convenience function for comprehensive validation
def validate_module_files(
    module_id: str,
    citation_content: Optional[str] = None,
    license_content: Optional[str] = None,
    omnibenchmark_content: Optional[str] = None,
    files_present: Optional[Dict[str, bool]] = None,
    warn_mode: bool = False,
    base_path: str = "",
) -> ValidationResult:
    """Validate all files for a module comprehensively."""
    ctx = create_validation_context(module_id, warn_mode, base_path)

    try:
        # First, check if CITATION.cff has a license field
        citation_has_license = False
        if citation_content is not None:
            try:
                citation_data = yaml.safe_load(citation_content)
                if isinstance(citation_data, dict):
                    citation_has_license = bool(citation_data.get("license"))
            except (yaml.YAMLError, AttributeError):
                pass

        # Validate file structure
        if files_present:
            structure_result = validate_file_structure(
                files_present, ctx, citation_has_license
            )
            ctx.merge_result(structure_result)

        # Validate CITATION.cff
        if citation_content is not None:
            citation_result, citation_data = validate_citation_cff_content(
                citation_content, ctx, "CITATION.cff"
            )
            ctx.merge_result(citation_result)

        # Validate LICENSE file
        if license_content is not None:
            license_result = validate_license_file_content(
                license_content, ctx, "LICENSE"
            )
            ctx.merge_result(license_result)

        # Validate omnibenchmark.yaml
        if omnibenchmark_content is not None:
            omnibenchmark_result, omnibenchmark_data = (
                validate_omnibenchmark_yaml_content(
                    omnibenchmark_content, ctx, "omnibenchmark.yaml"
                )
            )
            ctx.merge_result(omnibenchmark_result)

        # Validate license consistency
        if citation_content is not None and license_content is not None:
            try:
                citation_data = yaml.safe_load(citation_content)
                citation_license = (
                    citation_data.get("license")
                    if isinstance(citation_data, dict)
                    else None
                )
                consistency_result = validate_license_consistency(
                    citation_license, license_content, ctx
                )
                ctx.merge_result(consistency_result)
            except (yaml.YAMLError, AttributeError):
                pass  # Skip consistency check if citation is invalid

    except ValidationException as e:
        # In strict mode, return partial results with the exception issues
        ctx.result.add_issues(e.issues)

    return ctx.result
