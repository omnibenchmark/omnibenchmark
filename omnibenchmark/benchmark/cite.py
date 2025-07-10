"""Citation metadata extraction from benchmark modules."""

import logging
import yaml
import tempfile
import shutil
from pathlib import Path
from typing import Dict, Optional, Any, List

from omnibenchmark.benchmark.benchmark import Benchmark
from omnibenchmark.io.code import clone_module
from omnibenchmark.model.repo import get_repo_hash
from omnibenchmark.benchmark.validation_core import (
    ValidationResult,
    ValidationIssue,
    ValidationSeverity,
    ValidationException,
    create_validation_context,
    validate_citation_cff_content,
)

logger = logging.getLogger(__name__)


class CitationExtractionError(Exception):
    """Raised when citation extraction fails in strict mode."""

    def __init__(
        self, message: str, failed_modules: List[str], issues: List[ValidationIssue]
    ):
        super().__init__(message)
        self.failed_modules = failed_modules
        self.issues = issues


class ModuleCitationResult:
    """Result of citation extraction for a single module."""

    def __init__(self, module_id: str):
        self.module_id = module_id
        self.repository_url: Optional[str] = None
        self.commit_hash: Optional[str] = None
        self.citation_data: Optional[Dict[str, Any]] = None
        self.citation_file_found: bool = False
        self.local_repo_exists: bool = False
        self.validation_result: ValidationResult = ValidationResult()

    def add_issue(self, issue: ValidationIssue):
        """Add a validation issue."""
        self.validation_result.add_issue(issue)

    def is_valid(self) -> bool:
        """Return True if no errors."""
        return self.validation_result.is_valid()

    def has_warnings(self) -> bool:
        """Return True if warnings exist."""
        return self.validation_result.has_warnings()


def convert_errors_to_warnings(result: ModuleCitationResult) -> ModuleCitationResult:
    """Convert all errors to warnings for warn mode."""
    if not result.validation_result.errors:
        return result

    # Create new ValidationResult with all issues as warnings
    new_result = ValidationResult()
    for issue in result.validation_result.issues:
        warning_issue = ValidationIssue(
            issue_type=issue.issue_type,
            severity=ValidationSeverity.WARNING,
            path=issue.path,
            module_id=issue.module_id,
            message=issue.message,
            **issue.context,
        )
        new_result.add_issue(warning_issue)

    result.validation_result = new_result
    return result


def create_cite_issue(
    issue_type: str,
    module_id: str,
    message: str,
    severity: ValidationSeverity = ValidationSeverity.ERROR,
    **context,
) -> ValidationIssue:
    """Create a citation-related validation issue."""
    return ValidationIssue(
        issue_type=issue_type,
        severity=severity,
        path=".",
        module_id=module_id,
        message=message,
        **context,
    )


def extract_citation_metadata(
    benchmark: Benchmark, strict: bool = True, warn_mode: bool = False
) -> Dict[str, Optional[Dict[str, Any]]]:
    """Extract CITATION.cff metadata from all modules in a benchmark.

    Args:
        benchmark: The benchmark instance to extract citations from
        strict: If True, use strict mode (fail fast on first error)
        warn_mode: If True, convert errors to warnings instead of failing

    Returns:
        Dictionary mapping module_id to citation metadata

    Raises:
        RuntimeWarning: If no cloned repositories are found
        CitationExtractionError: If strict=True and errors found (unless warn_mode=True)
    """
    repos_base_dir = Path(".snakemake/repos")
    temp_directories = []  # Track temporary directories for cleanup

    modules = benchmark.get_converter().get_modules()
    module_results = {}
    found_any_repos = False
    all_issues = []

    try:
        for module_id, module in modules.items():
            result = _extract_single_module_citation(
                module_id, module, benchmark, repos_base_dir, temp_directories
            )

            if result.local_repo_exists:
                found_any_repos = True

            # In strict mode, fail fast on first error
            if strict and not warn_mode and not result.is_valid():
                all_issues.extend(result.validation_result.errors)
                raise CitationExtractionError(
                    f"Citation extraction failed for module: {module_id}",
                    [module_id],
                    result.validation_result.errors,
                )

            # Convert errors to warnings if in warn mode
            if warn_mode and not result.is_valid():
                result = convert_errors_to_warnings(result)
                # Log warnings
                for warning in result.validation_result.warnings:
                    logger.warning(warning.message)

            # Collect all errors for non-strict mode
            if not result.is_valid():
                all_issues.extend(result.validation_result.errors)

            module_results[module_id] = result

        # Check if we found any repositories. If there are no temporary directories, it might be because
        # we could not clone the repos from the given URLs.
        if not found_any_repos and not temp_directories:
            raise RuntimeWarning(
                "No local repositories found for any modules. Try running the benchmark first to clone the modules."
            )

        # In non-warn mode, check if we have any failures
        if not warn_mode and all_issues:
            failed_modules = [
                mid for mid, result in module_results.items() if not result.is_valid()
            ]
            raise CitationExtractionError(
                f"Missing valid citation metadata in {len(failed_modules)} modules",
                failed_modules,
                all_issues,
            )

        citation_metadata = {}
        for module_id, result in module_results.items():
            citation_metadata[module_id] = {
                "repository_url": result.repository_url,
                "commit_hash": result.commit_hash,
                "citation_data": result.citation_data,
                "citation_file_found": result.citation_file_found,
            }

        return citation_metadata

    finally:
        # Always cleanup temporary directories
        for temp_dir in temp_directories:
            try:
                if temp_dir.exists():
                    shutil.rmtree(temp_dir)
                    logger.debug(f"Cleaned up temporary directory: {temp_dir}")
            except Exception as e:
                logger.warning(f"Failed to cleanup temporary directory {temp_dir}: {e}")


def _extract_single_module_citation(
    module_id: str,
    module: Any,
    benchmark: Benchmark,
    repos_base_dir: Path,
    temp_directories: List[Path],
) -> ModuleCitationResult:
    """Extract citation metadata for a single module."""
    result = ModuleCitationResult(module_id)

    repo_info = benchmark.get_converter().get_module_repository(module)

    # Get repository URL and commit hash
    repo_url = getattr(repo_info, "url", None) or getattr(repo_info, "repository", None)
    commit_hash = (
        getattr(repo_info, "commit_hash", None)
        or getattr(repo_info, "commit", None)
        or getattr(repo_info, "version", None)
    )

    result.repository_url = repo_url
    result.commit_hash = commit_hash

    if not repo_url or not commit_hash:
        issue = create_cite_issue(
            "missing_repository_info",
            module_id,
            f"Missing repository info for module: {module_id}",
        )
        result.add_issue(issue)
        return result

    # Generate the expected local folder name
    folder_name = get_repo_hash(repo_url, commit_hash)
    local_repo_path = repos_base_dir / folder_name

    if not local_repo_path.exists():
        # Try to clone the repository to a temporary directory
        logger.info(
            f"Local repository not found for {module_id}, attempting to clone..."
        )
        temp_repo_path = _clone_to_temp(repo_url, commit_hash, module_id)

        if temp_repo_path is None:
            issue = create_cite_issue(
                "clone_failed",
                module_id,
                f"Failed to clone repository for {module_id}: {repo_url}",
                repo_url=repo_url,
            )
            result.add_issue(issue)
            return result

        # Use the temporary path instead and track for cleanup
        local_repo_path = temp_repo_path
        temp_directories.append(temp_repo_path.parent)

    result.local_repo_exists = True
    citation_file = local_repo_path / "CITATION.cff"

    if not citation_file.exists():
        issue = create_cite_issue(
            "citation_missing", module_id, "CITATION.cff file not found"
        )
        result.add_issue(issue)
        return result

    try:
        with open(citation_file, "r", encoding="utf-8") as f:
            citation_content = f.read()

        # Use structured validation
        # Create validation context for citation validation
        ctx = create_validation_context(module_id, warn_mode=True)
        try:
            citation_result, citation_data = validate_citation_cff_content(
                citation_content, ctx
            )
            # Merge citation validation results
            result.validation_result.add_issues(ctx.result.issues)
        except ValidationException as ve:
            # In citation validation, add the validation issues
            result.validation_result.add_issues(ve.issues)

        # Copy validation results
        # Citation validation results are already merged above

        if citation_data:
            result.citation_data = citation_data
            result.citation_file_found = True
            logger.info(f"Found valid CITATION.cff for module: {module_id}")

    except Exception as e:
        issue = create_cite_issue(
            "citation_read_error",
            module_id,
            f"Error reading CITATION.cff for {module_id}: {str(e)}",
            error=str(e),
        )
        result.add_issue(issue)

    return result


def convert_to_bibtex(citation_metadata: Dict[str, Optional[Dict[str, Any]]]) -> str:
    """Convert citation metadata to BibTeX format.

    Args:
        citation_metadata: Dictionary of citation metadata from extract_citation_metadata

    Returns:
        BibTeX formatted string
    """
    bibtex_entries = []

    for module_id, metadata in citation_metadata.items():
        if (
            not metadata
            or not metadata.get("citation_data")
            or not metadata.get("citation_file_found")
        ):
            continue

        citation_data = metadata["citation_data"]

        # Create BibTeX entry
        entry_type = "misc"  # Default type for software
        entry_id = module_id.replace("-", "_").replace(".", "_")

        bibtex_entry = f"@{entry_type}{{{entry_id},\n"

        # Add title
        if "title" in citation_data:
            bibtex_entry += f"  title = {{{citation_data['title']}}},\n"

        # Add authors
        if "authors" in citation_data:
            authors = _format_bibtex_authors(citation_data["authors"])
            if authors:
                bibtex_entry += f"  author = {{{authors}}},\n"

        # Add DOI
        if "doi" in citation_data:
            bibtex_entry += f"  doi = {{{citation_data['doi']}}},\n"

        # Add URL
        if "repository-code" in citation_data:
            bibtex_entry += f"  url = {{{citation_data['repository-code']}}},\n"
        elif metadata.get("repository_url"):
            bibtex_entry += f"  url = {{{metadata['repository_url']}}},\n"

        # Add year
        if "date-released" in citation_data:
            year = _extract_year_from_date(citation_data["date-released"])
            if year:
                bibtex_entry += f"  year = {{{year}}},\n"

        # Add version
        if "version" in citation_data:
            bibtex_entry += f"  version = {{{citation_data['version']}}},\n"

        # Add note about software
        bibtex_entry += "  note = {Software},\n"

        bibtex_entry = bibtex_entry.rstrip(",\n") + "\n}\n"
        bibtex_entries.append(bibtex_entry)

    return "\n".join(bibtex_entries) if bibtex_entries else "% No citation data found"


def format_output(
    citation_metadata: Dict[str, Optional[Dict[str, Any]]], format_type: str
) -> str:
    """Format citation metadata for output, avoiding YAML object serialization issues.

    Args:
        citation_metadata: Dictionary of citation metadata
        format_type: Output format ('json', 'yaml', 'bibtex')

    Returns:
        Formatted string
    """
    if format_type.lower() == "json":
        import json

        # Convert to simple dict to avoid serialization issues
        clean_data = _clean_for_serialization(citation_metadata)
        return json.dumps(clean_data, indent=2, default=str)
    elif format_type.lower() == "yaml":
        # Convert to simple dict to avoid YAML object serialization issues
        clean_data = _clean_for_serialization(citation_metadata)
        return yaml.dump(clean_data, default_flow_style=False)
    elif format_type.lower() == "bibtex":
        return convert_to_bibtex(citation_metadata)
    else:
        raise ValueError(f"Unsupported format: {format_type}")


def _clean_for_serialization(
    data: Dict[str, Optional[Dict[str, Any]]],
) -> Dict[str, Any]:
    """Clean data structure for serialization, removing complex objects."""
    clean_data = {}

    for module_id, metadata in data.items():
        if metadata is None:
            clean_data[str(module_id)] = None
            continue

        clean_metadata = {}
        for key, value in metadata.items():
            if key == "citation_data" and isinstance(value, dict):
                # Clean citation data recursively
                clean_metadata[key] = _clean_citation_data(value)
            else:
                clean_metadata[key] = str(value) if value is not None else None

        clean_data[str(module_id)] = clean_metadata

    return clean_data


def _clone_to_temp(repo_url: str, commit_hash: str, module_id: str) -> Optional[Path]:
    """Clone repository to temporary directory for citation extraction.

    Args:
        repo_url: Repository URL to clone
        commit_hash: Git commit hash or tag
        module_id: Module identifier for logging

    Returns:
        Path to temporary repository directory, or None if clone failed
    """
    try:
        # Create temporary directory
        temp_dir = Path(tempfile.mkdtemp(prefix=f"omnibenchmark_cite_{module_id}_"))
        logger.debug(f"Created temporary directory: {temp_dir}")

        # Clone the repository
        cloned_path = clone_module(temp_dir, repo_url, commit_hash)
        logger.info(
            f"Successfully cloned {repo_url}@{commit_hash} to temporary location"
        )

        return cloned_path

    except Exception as e:
        logger.warning(f"Failed to clone repository {repo_url}@{commit_hash}: {e}")
        # Clean up temp directory if it was created
        if "temp_dir" in locals() and temp_dir.exists():
            try:
                shutil.rmtree(temp_dir)
            except Exception as cleanup_error:
                logger.debug(
                    f"Failed to cleanup temp directory {temp_dir}: {cleanup_error}"
                )
        return None


def _clean_citation_data(citation_data: Dict[str, Any]) -> Dict[str, Any]:
    """Clean citation data for serialization."""
    clean_data = {}

    for key, value in citation_data.items():
        if isinstance(value, list):
            clean_data[key] = [_clean_citation_item(item) for item in value]
        elif isinstance(value, dict):
            clean_data[key] = _clean_citation_item(value)
        else:
            clean_data[key] = str(value) if value is not None else None

    return clean_data


def _clean_citation_item(item: Any) -> Any:
    """Clean individual citation items."""
    if isinstance(item, dict):
        return {k: str(v) if v is not None else None for k, v in item.items()}
    else:
        return str(item) if item is not None else None


def _format_bibtex_authors(authors: List[Dict]) -> str:
    """Format authors list for BibTeX."""
    formatted_authors = []

    for author in authors:
        if not isinstance(author, dict):
            continue

        family = author.get("family-names", "")
        given = author.get("given-names", "")

        if family and given:
            formatted_authors.append(f"{family}, {given}")
        elif family:
            formatted_authors.append(family)

    return " and ".join(formatted_authors)


def _extract_year_from_date(date_str) -> Optional[str]:
    """Extract year from date string."""
    try:
        return str(date_str).split("-")[0]
    except (AttributeError, IndexError):
        return None


def get_citation_summary(
    citation_metadata: Dict[str, Optional[Dict[str, Any]]],
) -> Dict[str, float]:
    """Get summary statistics of citation metadata.

    Args:
        citation_metadata: Dictionary of citation metadata

    Returns:
        Dictionary with summary statistics
    """
    total_modules = len(citation_metadata)
    found_citations = sum(
        1
        for m in citation_metadata.values()
        if m is not None and m.get("citation_file_found")
    )
    local_repos = sum(
        1
        for m in citation_metadata.values()
        if m is not None and m.get("local_repo_exists")
    )

    return {
        "total_modules": total_modules,
        "modules_with_citations": found_citations,
        "modules_with_local_repos": local_repos,
        "citation_coverage": found_citations / total_modules
        if total_modules > 0
        else 0.0,
    }


def cleanup_temp_repositories():
    """Clean up any temporary repositories created during citation extraction.

    This function removes temporary directories that may have been left behind
    due to interrupted execution or exceptions.
    """
    import glob

    temp_pattern = f"{tempfile.gettempdir()}/omnibenchmark_cite_*"

    for temp_path in glob.glob(temp_pattern):
        try:
            shutil.rmtree(temp_path)
            logger.debug(f"Cleaned up temporary directory: {temp_path}")
        except Exception as e:
            logger.warning(f"Failed to cleanup temporary directory {temp_path}: {e}")
