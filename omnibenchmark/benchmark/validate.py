"""Comprehensive validation functionality for benchmark modules."""

import logging
import yaml
from pathlib import Path
from typing import Dict, Optional, Any, List

from omnibenchmark.benchmark.benchmark import BenchmarkExecution
from omnibenchmark.benchmark.repository_utils import (
    RepositoryManager,
    get_module_repository_info,
    resolve_module_repository,
)
from omnibenchmark.model.module import ModuleMetadata
import spdx_license_list

logger = logging.getLogger(__name__)


class ValidationError(Exception):
    """Raised when validation fails in strict mode."""

    def __init__(self, message: str, failed_modules: List[str]):
        super().__init__(message)
        self.failed_modules = failed_modules


class ValidationResult:
    """Container for comprehensive validation results."""

    def __init__(self, module_id: str):
        self.module_id = module_id
        self.repository_url: Optional[str] = None
        self.commit_hash: Optional[str] = None
        self.local_repo_exists: bool = False

        # CITATION.cff validation
        self.citation_file_exists: bool = False
        self.citation_file_valid: bool = False
        self.citation_has_license: bool = False
        self.citation_has_authors: bool = False
        self.citation_license: Optional[str] = None
        self.citation_data: Optional[Dict[str, Any]] = None

        # License validation
        self.license_file_exists: bool = False
        self.license_file_content: Optional[str] = None
        self.license_matches_citation: bool = False
        self.detected_license: Optional[str] = None

        # Module metadata validation
        self.omnibenchmark_yaml_exists: bool = False
        self.omnibenchmark_yaml_valid: bool = False
        self.module_metadata: Optional[ModuleMetadata] = None

        self.validation_errors: List[str] = []
        self.validation_warnings: List[str] = []

    def is_valid(self) -> bool:
        """Check if module passed all validations."""
        return (
            self.repository_url is not None
            and self.commit_hash is not None
            and self.local_repo_exists
            and self.citation_file_exists
            and self.citation_file_valid
            and self.citation_has_license
            and self.citation_has_authors
            and len(self.validation_errors) == 0
        )

    def add_error(self, error: str):
        """Add a validation error."""
        self.validation_errors.append(error)

    def add_warning(self, warning: str):
        """Add a validation warning."""
        self.validation_warnings.append(warning)


def validate_modules(
    benchmark: BenchmarkExecution, strict: bool = False
) -> Dict[str, ValidationResult]:
    """Comprehensively validate all benchmark modules.

    Validates:
    - CITATION.cff existence and format
    - Required CITATION.cff fields (license, authors)
    - LICENSE file existence and consistency with CITATION.cff
    - omnibenchmark.yaml file existence and validity

    Args:
        benchmark: The benchmark instance to validate
        strict: If True, raise errors for ANY module with validation issues

    Returns:
        Dictionary mapping module_id to ValidationResult

    Raises:
        RuntimeWarning: If no cloned repositories are found
        ValidationError: If strict=True and ANY module has validation issues
    """
    modules = benchmark.get_model().get_modules()
    validation_results = {}
    found_any_repos = False
    failed_modules = []

    with RepositoryManager(prefix="omnibenchmark_validate") as repo_manager:
        for module_id, module in modules.items():
            result = ValidationResult(module_id)
            validation_results[module_id] = result

            # Get repository information
            repo_url, commit_hash = get_module_repository_info(benchmark, module)

            result.repository_url = repo_url
            result.commit_hash = commit_hash

            if not repo_url or not commit_hash:
                result.add_error("Missing repository info (URL or commit hash)")
                if strict:
                    failed_modules.append(module_id)
                else:
                    logger.warning(f"Missing repository info for module: {module_id}")
                continue

            # Resolve repository path (local or clone to temp)
            local_repo_path = resolve_module_repository(
                benchmark, module, module_id, repo_manager
            )

            if local_repo_path is None:
                result.add_error("Failed to access repository")
                if strict:
                    failed_modules.append(module_id)
                else:
                    logger.warning(
                        f"Failed to access repository for {module_id}: {repo_url}"
                    )
                continue

            found_any_repos = True
            result.local_repo_exists = True

            # Comprehensive validation of the module
            _validate_module_comprehensive(local_repo_path, result)

            if result.validation_errors:
                if strict:
                    failed_modules.append(module_id)
                else:
                    for error in result.validation_errors:
                        logger.warning(f"Validation error in {module_id}: {error}")

            if result.is_valid():
                logger.info(f"Successfully validated module: {module_id}")

        # Check if we found any repositories
        if not found_any_repos:
            raise RuntimeWarning(
                "No local repositories found for any modules. Try running the benchmark first to clone the modules."
            )

        if strict and failed_modules:
            raise ValidationError(
                f"Validation failed for {len(failed_modules)} modules", failed_modules
            )

        return validation_results


def _validate_module_comprehensive(repo_path: Path, result: ValidationResult):
    """Perform comprehensive validation of a module repository.

    Args:
        repo_path: Path to the cloned repository
        result: ValidationResult to populate
    """
    # 1. Validate CITATION.cff
    _validate_citation_cff(repo_path, result)

    # 2. Validate LICENSE file
    _validate_license_file(repo_path, result)

    # 3. Check license consistency
    _validate_license_consistency(result)

    # 4. Validate omnibenchmark.yaml
    _validate_omnibenchmark_yaml(repo_path, result)


def _validate_citation_cff(repo_path: Path, result: ValidationResult):
    """Validate CITATION.cff file."""
    citation_file = repo_path / "CITATION.cff"

    if not citation_file.exists():
        result.add_error(
            f"CITATION.cff file not found in repository at {repo_path} (module: {result.module_id})"
        )
        return

    result.citation_file_exists = True

    try:
        with open(citation_file, "r", encoding="utf-8") as f:
            citation_data = yaml.safe_load(f)
        result.citation_data = citation_data
    except yaml.YAMLError as e:
        result.add_error(f"Invalid YAML syntax in CITATION.cff: {e}")
        return
    except Exception as e:
        result.add_error(f"Failed to read CITATION.cff: {e}")
        return

    if not isinstance(citation_data, dict):
        result.add_error("CITATION.cff must contain a YAML object")
        return

    # Validate required fields
    required_fields = ["cff-version", "message", "title"]
    for field in required_fields:
        if field not in citation_data:
            result.add_error(f"Missing required field in CITATION.cff: {field}")

    # Check for authors (REQUIRED)
    if "authors" not in citation_data or not citation_data["authors"]:
        result.add_error("CITATION.cff must have authors field")
    else:
        result.citation_has_authors = True
        authors = citation_data["authors"]
        if not isinstance(authors, list):
            result.add_error("Authors field must be a list")
        else:
            for i, author in enumerate(authors):
                if not isinstance(author, dict):
                    result.add_error(f"Author {i + 1}: must be an object")
                    continue
                if "family-names" not in author:
                    result.add_error(
                        f"Author {i + 1}: missing required field 'family-names'"
                    )

    # Check for license (REQUIRED)
    if "license" not in citation_data or not citation_data["license"]:
        result.add_error("CITATION.cff must have license field")
    else:
        result.citation_has_license = True
        result.citation_license = str(citation_data["license"])

        # Validate license using SPDX
        if result.citation_license not in spdx_license_list.LICENSES.keys():
            result.add_error(
                f"License '{result.citation_license}' not in recognized SPDX licenses"
            )

    # Validate CFF version
    if "cff-version" in citation_data:
        cff_version = str(citation_data["cff-version"])
        if not cff_version.startswith("1."):
            result.add_error(f"Unsupported CFF version: {cff_version}")

    if not result.validation_errors:
        result.citation_file_valid = True


def _validate_license_file(repo_path: Path, result: ValidationResult):
    """Validate LICENSE file existence and content."""
    # Check for common LICENSE file names
    license_files = ["LICENSE", "LICENSE.txt", "LICENSE.md", "COPYING", "COPYING.txt"]
    license_file_path = None

    for filename in license_files:
        potential_path = repo_path / filename
        if potential_path.exists():
            license_file_path = potential_path
            break

    if license_file_path:
        result.license_file_exists = True
        try:
            with open(license_file_path, "r", encoding="utf-8") as f:
                result.license_file_content = f.read()
        except Exception as e:
            result.add_warning(f"Could not read LICENSE file: {e}")
    else:
        # If no license in CITATION.cff, missing LICENSE file is an error
        if not result.citation_has_license:
            result.add_error("No LICENSE file found and no license in CITATION.cff")
        else:
            result.add_warning("No LICENSE file found in repository")


def _validate_license_consistency(result: ValidationResult):
    """Validate that CITATION.cff license matches LICENSE file."""
    if not result.citation_has_license:
        return  # Already flagged as error

    if not result.license_file_exists:
        return  # Already handled in _validate_license_file

    if not result.license_file_content or not result.citation_license:
        return

    # Simple license detection in LICENSE file content
    detected_license = _detect_license_from_content(result.license_file_content)
    result.detected_license = detected_license

    if detected_license and detected_license != result.citation_license:
        result.add_error(
            f"License mismatch - CITATION.cff: '{result.citation_license}' vs LICENSE file: '{detected_license}'"
        )
    elif detected_license:
        result.license_matches_citation = True
    else:
        result.add_warning(
            "Could not automatically detect license type from LICENSE file"
        )


def _validate_omnibenchmark_yaml(repo_path: Path, result: ValidationResult):
    """Validate omnibenchmark.yaml file or legacy config.cfg."""
    omnibenchmark_file = repo_path / "omnibenchmark.yaml"
    legacy_config_file = repo_path / "config.cfg"

    if not omnibenchmark_file.exists():
        # Check for legacy config.cfg
        if legacy_config_file.exists():
            result.omnibenchmark_yaml_exists = True
            result.omnibenchmark_yaml_valid = True
            result.add_warning(
                "Using legacy config.cfg file. Please migrate to omnibenchmark.yaml"
            )
            return
        else:
            result.add_warning("omnibenchmark.yaml file not found")
            return

    result.omnibenchmark_yaml_exists = True

    try:
        with open(omnibenchmark_file, "r", encoding="utf-8") as f:
            omnibenchmark_data = yaml.safe_load(f)

        # Basic validation - just check it's valid YAML and is a dict
        if isinstance(omnibenchmark_data, dict):
            result.omnibenchmark_yaml_valid = True
        else:
            result.add_warning("omnibenchmark.yaml should contain a YAML object")

    except yaml.YAMLError as e:
        result.add_warning(f"Invalid YAML syntax in omnibenchmark.yaml: {e}")
    except Exception as e:
        result.add_warning(f"Failed to read omnibenchmark.yaml: {e}")


def _detect_license_from_content(content: str) -> Optional[str]:
    """Detect license type from LICENSE file content."""
    content_lower = content.lower()

    # Simple pattern matching for common licenses
    # TODO: extend with more licenses
    license_patterns = {
        "MIT": ["mit license", "permission is hereby granted, free of charge"],
        "Apache-2.0": ["apache license", "version 2.0"],
        "GPL-3.0": ["gnu general public license", "version 3"],
        "GPL-2.0": ["gnu general public license", "version 2"],
        "BSD-3-Clause": [
            "bsd 3-clause",
            "redistribution and use in source and binary forms",
        ],
        "BSD-2-Clause": [
            "bsd 2-clause",
            "redistribution and use in source and binary forms",
        ],
    }

    for license_id, patterns in license_patterns.items():
        if all(pattern in content_lower for pattern in patterns):
            return license_id

    return None


def get_validation_summary(
    validation_results: Dict[str, ValidationResult],
) -> Dict[str, Any]:
    """Get summary statistics of comprehensive module validation.

    Args:
        validation_results: Dictionary of validation results

    Returns:
        Dictionary with summary statistics
    """
    total_modules = len(validation_results)
    valid_modules = sum(
        1 for result in validation_results.values() if result.is_valid()
    )
    modules_with_repos = sum(
        1 for result in validation_results.values() if result.local_repo_exists
    )
    modules_with_citation_files = sum(
        1 for result in validation_results.values() if result.citation_file_exists
    )
    modules_with_valid_citation = sum(
        1 for result in validation_results.values() if result.citation_file_valid
    )
    modules_with_license = sum(
        1 for result in validation_results.values() if result.citation_has_license
    )
    modules_with_authors = sum(
        1 for result in validation_results.values() if result.citation_has_authors
    )
    modules_with_license_file = sum(
        1 for result in validation_results.values() if result.license_file_exists
    )
    modules_with_omnibenchmark_yaml = sum(
        1 for result in validation_results.values() if result.omnibenchmark_yaml_exists
    )

    return {
        "total_modules": total_modules,
        "valid_modules": valid_modules,
        "modules_with_repos": modules_with_repos,
        "modules_with_citation_files": modules_with_citation_files,
        "modules_with_valid_citation": modules_with_valid_citation,
        "modules_with_license": modules_with_license,
        "modules_with_authors": modules_with_authors,
        "modules_with_license_file": modules_with_license_file,
        "modules_with_omnibenchmark_yaml": modules_with_omnibenchmark_yaml,
        "validation_success_rate": valid_modules / total_modules
        if total_modules > 0
        else 0.0,
    }


def format_validation_results(
    validation_results: Dict[str, ValidationResult], format_type: str = "summary"
) -> str:
    """Format comprehensive validation results for output.

    Args:
        validation_results: Dictionary of validation results
        format_type: Output format ('summary', 'detailed', 'json')

    Returns:
        Formatted string
    """
    if format_type == "json":
        import json

        data = {}
        for module_id, result in validation_results.items():
            data[module_id] = {
                "valid": result.is_valid(),
                "repository_url": result.repository_url,
                "commit_hash": result.commit_hash,
                "local_repo_exists": result.local_repo_exists,
                "citation_file_exists": result.citation_file_exists,
                "citation_file_valid": result.citation_file_valid,
                "citation_has_license": result.citation_has_license,
                "citation_has_authors": result.citation_has_authors,
                "license_file_exists": result.license_file_exists,
                "license_matches_citation": result.license_matches_citation,
                "omnibenchmark_yaml_exists": result.omnibenchmark_yaml_exists,
                "omnibenchmark_yaml_valid": result.omnibenchmark_yaml_valid,
                "errors": result.validation_errors,
                "warnings": result.validation_warnings,
            }
        return json.dumps(data, indent=2)

    elif format_type == "detailed":
        output = []
        for module_id, result in validation_results.items():
            output.append(f"\n=== Module: {module_id} ===")
            output.append(f"Overall Valid: {'✅' if result.is_valid() else '❌'}")
            output.append(f"Repository URL: {result.repository_url or 'N/A'}")
            output.append(f"Commit Hash: {result.commit_hash or 'N/A'}")

            output.append("\nCITATION.cff:")
            output.append(f"  Exists: {'✅' if result.citation_file_exists else '❌'}")
            output.append(f"  Valid: {'✅' if result.citation_file_valid else '❌'}")
            output.append(
                f"  Has License: {'✅' if result.citation_has_license else '❌'}"
            )
            output.append(
                f"  Has Authors: {'✅' if result.citation_has_authors else '❌'}"
            )

            output.append("\nLicense:")
            output.append(
                f"  LICENSE file exists: {'✅' if result.license_file_exists else '❌'}"
            )

            # Only show license matching if CITATION.cff exists
            if result.citation_file_exists:
                output.append(
                    f"  License matches CITATION.cff: {'✅' if result.license_matches_citation else '❌'}"
                )
                if result.citation_license:
                    output.append(f"  CITATION.cff license: {result.citation_license}")
                if result.detected_license:
                    output.append(f"  Detected LICENSE file: {result.detected_license}")
                if (
                    result.citation_license
                    and result.detected_license
                    and result.citation_license != result.detected_license
                ):
                    output.append("  ⚠️ License mismatch detected!")
            else:
                if result.detected_license:
                    output.append(f"  Detected LICENSE file: {result.detected_license}")

            output.append("\nOmnibenchmark:")
            output.append(
                f"  omnibenchmark.yaml exists: {'✅' if result.omnibenchmark_yaml_exists else '❌'}"
            )
            output.append(
                f"  omnibenchmark.yaml valid: {'✅' if result.omnibenchmark_yaml_valid else '❌'}"
            )

            if result.validation_errors:
                output.append("\nErrors:")
                for error in result.validation_errors:
                    output.append(f"  ❌ {error}")

            if result.validation_warnings:
                output.append("\nWarnings:")
                for warning in result.validation_warnings:
                    output.append(f"  ⚠ {warning}")

        return "\n".join(output)

    else:  # summary
        summary = get_validation_summary(validation_results)
        valid = summary["valid_modules"]
        total = summary["total_modules"]
        success_rate = summary["validation_success_rate"] * 100

        output = [
            "Comprehensive Module Validation Summary:",
            f"  Total modules: {total}",
            f"  Fully valid modules: {valid}",
            f"  Success rate: {success_rate:.1f}%",
            "",
            "Component Status:",
            f"  Modules with repositories: {summary['modules_with_repos']}",
            f"  Modules with CITATION.cff: {summary['modules_with_citation_files']}",
            f"  Modules with valid CITATION.cff: {summary['modules_with_valid_citation']}",
            f"  Modules with license in CITATION.cff: {summary['modules_with_license']}",
            f"  Modules with authors in CITATION.cff: {summary['modules_with_authors']}",
            f"  Modules with LICENSE file: {summary['modules_with_license_file']}",
            f"  Modules with omnibenchmark.yaml: {summary['modules_with_omnibenchmark_yaml']}",
        ]

        # Add failed modules
        failed_modules = [
            module_id
            for module_id, result in validation_results.items()
            if not result.is_valid()
        ]
        if failed_modules:
            output.append(f"\nFailed modules: {', '.join(failed_modules)}")

        return "\n".join(output)
