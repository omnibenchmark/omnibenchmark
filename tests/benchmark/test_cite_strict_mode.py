import pytest

from omnibenchmark.benchmark.cite import (
    CitationExtractionError,
    ModuleCitationResult,
    convert_errors_to_warnings,
    create_cite_issue,
)
from omnibenchmark.benchmark.metadata import (
    ValidationSeverity,
)


@pytest.mark.short
def test_module_citation_result_init():
    """Test ModuleCitationResult initialization."""
    result = ModuleCitationResult("test_module")

    assert result.module_id == "test_module"
    assert result.repository_url is None
    assert result.commit_hash is None
    assert result.citation_data is None
    assert result.citation_file_found is False
    assert result.local_repo_exists is False
    assert result.is_valid() is True
    assert result.has_warnings() is False


@pytest.mark.short
def test_module_citation_result_add_issue():
    """Test adding issues to ModuleCitationResult."""
    result = ModuleCitationResult("test_module")
    issue = create_cite_issue(
        "missing_repository_info", "test_module", "Missing repository info"
    )

    result.add_issue(issue)

    assert not result.is_valid()
    assert len(result.validation_result.errors) == 1
    assert result.validation_result.errors[0].issue_type == issue.issue_type


@pytest.mark.short
def test_module_citation_result_add_warning():
    """Test adding warnings to ModuleCitationResult."""
    result = ModuleCitationResult("test_module")
    warning = create_cite_issue(
        "missing_repository_info",
        "test_module",
        "Missing repository info",
        ValidationSeverity.WARNING,
    )

    result.add_issue(warning)

    assert result.is_valid()  # Warnings don't affect validity
    assert result.has_warnings()
    assert len(result.validation_result.warnings) == 1
    assert result.validation_result.warnings[0].issue_type == warning.issue_type


@pytest.mark.short
def test_create_cite_issue():
    """Test creating citation-related validation issues."""
    error_issue = create_cite_issue(
        "missing_repository_info", "test_module", "Missing repository info"
    )
    warning_issue = create_cite_issue(
        "missing_repository_info",
        "test_module",
        "Missing repository info",
        ValidationSeverity.WARNING,
    )

    assert error_issue.issue_type == "missing_repository_info"
    assert error_issue.module_id == "test_module"
    assert error_issue.severity == ValidationSeverity.ERROR
    assert error_issue.message == "Missing repository info"

    assert warning_issue.severity == ValidationSeverity.WARNING


@pytest.mark.short
def test_convert_errors_to_warnings():
    """Test converting all errors in a result to warnings."""
    result = ModuleCitationResult("test_module")

    # Add multiple errors
    result.add_issue(
        create_cite_issue(
            "missing_repository_info", "test_module", "Missing repository info"
        )
    )
    result.add_issue(
        create_cite_issue(
            "clone_failed",
            "test_module",
            "Clone failed",
            repo_url="https://github.com/test/repo.git",
        )
    )

    assert not result.is_valid()
    assert len(result.validation_result.errors) == 2
    assert len(result.validation_result.warnings) == 0

    # Convert errors to warnings
    converted_result = convert_errors_to_warnings(result)

    assert converted_result.is_valid()
    assert len(converted_result.validation_result.errors) == 0
    assert len(converted_result.validation_result.warnings) == 2

    # Check that all issues are now warnings
    for issue in converted_result.validation_result.warnings:
        assert issue.severity == ValidationSeverity.WARNING


@pytest.mark.short
def test_convert_errors_to_warnings_no_errors():
    """Test converting errors to warnings when no errors exist."""
    result = ModuleCitationResult("test_module")

    assert result.is_valid()

    # Convert should be no-op
    converted_result = convert_errors_to_warnings(result)

    assert converted_result.is_valid()
    assert len(converted_result.validation_result.errors) == 0
    assert len(converted_result.validation_result.warnings) == 0


@pytest.mark.short
def test_citation_extraction_error_attributes():
    """Test CitationExtractionError contains proper attributes."""
    failed_modules = ["module1", "module2"]
    issues = [
        create_cite_issue(
            "missing_repository_info", "module1", "Missing repository info"
        ),
        create_cite_issue(
            "clone_failed",
            "module2",
            "Clone failed",
            repo_url="https://github.com/test/repo.git",
        ),
    ]

    exc = CitationExtractionError("Test message", failed_modules, issues)

    assert str(exc) == "Test message"
    assert exc.failed_modules == failed_modules
    assert exc.issues == issues
    assert len(exc.issues) == 2


@pytest.mark.short
def test_issue_consistency():
    """Test that error and warning issues have consistent messages."""
    module_id = "test_module"
    repo_url = "https://github.com/test/repo.git"
    error_msg = "Test error"

    # Test MissingRepositoryInfo
    error = create_cite_issue(
        "missing_repository_info", module_id, "Missing repository info"
    )
    warning = create_cite_issue(
        "missing_repository_info",
        module_id,
        "Missing repository info",
        ValidationSeverity.WARNING,
    )
    assert error.message == warning.message

    # Test CloneFailed
    error = create_cite_issue(
        "clone_failed", module_id, "Clone failed", repo_url=repo_url
    )
    warning = create_cite_issue(
        "clone_failed",
        module_id,
        "Clone failed",
        ValidationSeverity.WARNING,
        repo_url=repo_url,
    )
    assert error.message == warning.message

    # Test CitationRead
    error = create_cite_issue(
        "citation_read_error", module_id, error_msg, error=error_msg
    )
    warning = create_cite_issue(
        "citation_read_error",
        module_id,
        error_msg,
        ValidationSeverity.WARNING,
        error=error_msg,
    )
    assert error.message == warning.message
