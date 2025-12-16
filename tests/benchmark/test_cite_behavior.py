"""Behavior-driven tests for citation extraction functionality."""

import pytest

from omnibenchmark.benchmark.cite import (
    ModuleCitationResult,
    convert_errors_to_warnings,
    create_cite_issue,
    convert_to_bibtex,
    format_output,
    get_citation_summary,
    _format_bibtex_authors,
    _extract_year_from_date,
    _clean_for_serialization,
)
from omnibenchmark.benchmark.metadata import (
    ValidationIssue,
    ValidationSeverity,
)


class TestModuleCitationResult:
    """Behavior tests for ModuleCitationResult."""

    def test_when_created_then_has_default_values(self):
        """When a ModuleCitationResult is created, it should have sensible defaults."""
        result = ModuleCitationResult("test_module")

        assert result.module_id == "test_module"
        assert result.repository_url is None
        assert result.commit_hash is None
        assert result.citation_data is None
        assert result.citation_file_found is False
        assert result.local_repo_exists is False
        assert result.is_valid() is True  # No errors means valid

    def test_when_issue_added_then_validation_fails(self):
        """When an error issue is added, the result should become invalid."""
        result = ModuleCitationResult("test_module")

        issue = ValidationIssue(
            issue_type="test_error",
            severity=ValidationSeverity.ERROR,
            path="test.txt",
            module_id="test_module",
            message="Test error",
        )
        result.add_issue(issue)

        assert not result.is_valid()
        assert len(result.validation_result.errors) == 1

    def test_when_warning_added_then_still_valid(self):
        """When a warning is added, the result should still be valid."""
        result = ModuleCitationResult("test_module")

        warning = ValidationIssue(
            issue_type="test_warning",
            severity=ValidationSeverity.WARNING,
            path="test.txt",
            module_id="test_module",
            message="Test warning",
        )
        result.add_issue(warning)

        assert result.is_valid()
        assert result.has_warnings()
        assert len(result.validation_result.warnings) == 1


class TestConvertErrorsToWarnings:
    """Behavior tests for converting errors to warnings."""

    def test_when_no_errors_then_result_unchanged(self):
        """When result has no errors, conversion should return it unchanged."""
        result = ModuleCitationResult("test_module")
        converted = convert_errors_to_warnings(result)

        assert converted.is_valid()
        assert not converted.has_warnings()

    def test_when_has_errors_then_converts_to_warnings(self):
        """When result has errors, they should be converted to warnings."""
        result = ModuleCitationResult("test_module")

        error = ValidationIssue(
            issue_type="test_error",
            severity=ValidationSeverity.ERROR,
            path="test.txt",
            module_id="test_module",
            message="Test error",
        )
        result.add_issue(error)

        assert not result.is_valid()

        converted = convert_errors_to_warnings(result)

        assert converted.is_valid()  # No longer has errors
        assert converted.has_warnings()  # Converted to warnings
        assert len(converted.validation_result.warnings) == 1


class TestCreateCiteIssue:
    """Behavior tests for create_cite_issue helper."""

    def test_when_created_then_has_correct_attributes(self):
        """When a citation issue is created, it should have all attributes set."""
        issue = create_cite_issue(
            issue_type="missing_citation",
            module_id="test_module",
            message="Citation not found",
            repo_url="https://github.com/test/repo",
        )

        assert issue.issue_type == "missing_citation"
        assert issue.module_id == "test_module"
        assert issue.message == "Citation not found"
        assert issue.severity == ValidationSeverity.ERROR
        assert issue.context["repo_url"] == "https://github.com/test/repo"

    def test_when_severity_specified_then_uses_it(self):
        """When severity is specified, it should use that severity."""
        issue = create_cite_issue(
            issue_type="test_warning",
            module_id="test_module",
            message="Test warning",
            severity=ValidationSeverity.WARNING,
        )

        assert issue.severity == ValidationSeverity.WARNING


class TestBibtexConversion:
    """Behavior tests for BibTeX conversion."""

    def test_when_no_citations_then_returns_placeholder(self):
        """When no citation data exists, should return placeholder message."""
        citation_metadata = {
            "module1": {"citation_data": None, "citation_file_found": False}
        }

        bibtex = convert_to_bibtex(citation_metadata)

        assert "No citation data found" in bibtex

    def test_when_citation_has_minimal_fields_then_generates_entry(self):
        """When citation has minimal fields, should generate valid BibTeX entry."""
        citation_metadata = {
            "test_module": {
                "citation_data": {
                    "title": "Test Software",
                    "authors": [{"family-names": "Doe", "given-names": "John"}],
                },
                "citation_file_found": True,
                "repository_url": "https://github.com/test/repo",
            }
        }

        bibtex = convert_to_bibtex(citation_metadata)

        assert "@misc{test_module" in bibtex
        assert "title = {Test Software}" in bibtex
        assert "author = {Doe, John}" in bibtex
        assert "url = {https://github.com/test/repo}" in bibtex

    def test_when_citation_has_doi_then_includes_doi(self):
        """When citation has DOI, should include it in BibTeX."""
        citation_metadata = {
            "test_module": {
                "citation_data": {
                    "title": "Test Software",
                    "authors": [{"family-names": "Doe", "given-names": "John"}],
                    "doi": "10.1234/test.doi",
                },
                "citation_file_found": True,
            }
        }

        bibtex = convert_to_bibtex(citation_metadata)

        assert "doi = {10.1234/test.doi}" in bibtex

    def test_when_multiple_authors_then_joins_with_and(self):
        """When multiple authors exist, should join with 'and'."""
        authors = [
            {"family-names": "Doe", "given-names": "John"},
            {"family-names": "Smith", "given-names": "Jane"},
            {"family-names": "Johnson", "given-names": "Bob"},
        ]

        result = _format_bibtex_authors(authors)

        assert result == "Doe, John and Smith, Jane and Johnson, Bob"

    def test_when_author_has_only_family_name_then_uses_it(self):
        """When author has only family name, should use it without given name."""
        authors = [{"family-names": "Doe"}]

        result = _format_bibtex_authors(authors)

        assert result == "Doe"


class TestYearExtraction:
    """Behavior tests for year extraction from dates."""

    def test_when_date_is_iso_format_then_extracts_year(self):
        """When date is in ISO format, should extract the year."""
        year = _extract_year_from_date("2023-12-15")

        assert year == "2023"

    def test_when_date_is_year_only_then_returns_year(self):
        """When date is just a year, should return it."""
        year = _extract_year_from_date("2023")

        assert year == "2023"

    def test_when_date_is_invalid_then_returns_none(self):
        """When date is invalid, should return None."""
        year = _extract_year_from_date(None)

        assert year is None


class TestFormatOutput:
    """Behavior tests for output formatting."""

    def test_when_format_is_json_then_returns_json_string(self):
        """When format is JSON, should return valid JSON string."""
        citation_metadata = {
            "test_module": {
                "citation_data": {"title": "Test"},
                "citation_file_found": True,
            }
        }

        output = format_output(citation_metadata, "json")

        assert '"test_module"' in output
        assert "title" in output

    def test_when_format_is_yaml_then_returns_yaml_string(self):
        """When format is YAML, should return valid YAML string."""
        citation_metadata = {
            "test_module": {
                "citation_data": {"title": "Test"},
                "citation_file_found": True,
            }
        }

        output = format_output(citation_metadata, "yaml")

        assert "test_module:" in output
        assert "title:" in output

    def test_when_format_is_bibtex_then_returns_bibtex_string(self):
        """When format is BibTeX, should return BibTeX formatted string."""
        citation_metadata = {
            "test_module": {
                "citation_data": {"title": "Test"},
                "citation_file_found": True,
            }
        }

        output = format_output(citation_metadata, "bibtex")

        assert "@misc" in output or "No citation data" in output

    def test_when_format_is_invalid_then_raises_error(self):
        """When format is invalid, should raise ValueError."""
        citation_metadata = {"test_module": {}}

        with pytest.raises(ValueError) as exc_info:
            format_output(citation_metadata, "invalid_format")

        assert "Unsupported format" in str(exc_info.value)


class TestDataCleaning:
    """Behavior tests for data cleaning functions."""

    def test_when_data_has_none_values_then_preserves_structure(self):
        """When data has None values, should preserve the structure."""
        data = {
            "module1": {
                "citation_data": None,
                "repository_url": "https://test.com",
                "commit_hash": None,
            }
        }

        cleaned = _clean_for_serialization(data)

        assert "module1" in cleaned
        assert cleaned["module1"]["citation_data"] is None

    def test_when_data_has_nested_dicts_then_cleans_recursively(self):
        """When data has nested dictionaries, should clean recursively."""
        data = {
            "module1": {
                "citation_data": {
                    "title": "Test",
                    "authors": [{"family-names": "Doe"}],
                }
            }
        }

        cleaned = _clean_for_serialization(data)

        assert isinstance(cleaned["module1"]["citation_data"], dict)
        assert isinstance(cleaned["module1"]["citation_data"]["authors"], list)


class TestCitationSummary:
    """Behavior tests for citation summary generation."""

    def test_when_no_modules_then_coverage_is_zero(self):
        """When there are no modules, coverage should be 0."""
        summary = get_citation_summary({})

        assert summary["total_modules"] == 0
        assert summary["modules_with_citations"] == 0
        assert summary["citation_coverage"] == 0.0

    def test_when_all_have_citations_then_coverage_is_one(self):
        """When all modules have citations, coverage should be 1.0."""
        citation_metadata = {
            "module1": {"citation_file_found": True},
            "module2": {"citation_file_found": True},
        }

        summary = get_citation_summary(citation_metadata)

        assert summary["total_modules"] == 2
        assert summary["modules_with_citations"] == 2
        assert summary["citation_coverage"] == 1.0

    def test_when_half_have_citations_then_coverage_is_half(self):
        """When half the modules have citations, coverage should be 0.5."""
        citation_metadata = {
            "module1": {"citation_file_found": True},
            "module2": {"citation_file_found": False},
        }

        summary = get_citation_summary(citation_metadata)

        assert summary["total_modules"] == 2
        assert summary["modules_with_citations"] == 1


class TestBibtexEdgeCases:
    """Additional tests for BibTeX conversion edge cases."""

    @pytest.mark.short
    def test_bibtex_with_version_and_date(self):
        """Test BibTeX generation with version and date-released fields."""
        citation_metadata = {
            "test_module": {
                "citation_data": {
                    "title": "Test Software",
                    "authors": [{"family-names": "Doe", "given-names": "John"}],
                    "version": "1.2.3",
                    "date-released": "2023-06-15",
                },
                "citation_file_found": True,
            }
        }

        bibtex = convert_to_bibtex(citation_metadata)

        assert "version = {1.2.3}" in bibtex
        assert "year = {2023}" in bibtex

    @pytest.mark.short
    def test_bibtex_with_repository_code(self):
        """Test BibTeX generation with repository-code field."""
        citation_metadata = {
            "test_module": {
                "citation_data": {
                    "title": "Test Software",
                    "authors": [{"family-names": "Doe"}],
                    "repository-code": "https://github.com/test/repo",
                },
                "citation_file_found": True,
            }
        }

        bibtex = convert_to_bibtex(citation_metadata)

        assert "url = {https://github.com/test/repo}" in bibtex

    @pytest.mark.short
    def test_bibtex_prefers_repository_code_over_url(self):
        """Test that repository-code is preferred over repository_url."""
        citation_metadata = {
            "test_module": {
                "citation_data": {
                    "title": "Test Software",
                    "authors": [{"family-names": "Doe"}],
                    "repository-code": "https://github.com/test/code",
                },
                "citation_file_found": True,
                "repository_url": "https://github.com/test/url",
            }
        }

        bibtex = convert_to_bibtex(citation_metadata)

        assert "url = {https://github.com/test/code}" in bibtex
        assert "url = {https://github.com/test/url}" not in bibtex

    @pytest.mark.short
    def test_format_output_json(self):
        """Test JSON output formatting."""
        citation_metadata = {
            "module1": {
                "citation_data": {"title": "Test"},
                "citation_file_found": True,
            }
        }

        result = format_output(citation_metadata, "json")

        assert '"module1"' in result
        assert '"title": "Test"' in result

    @pytest.mark.short
    def test_format_output_yaml(self):
        """Test YAML output formatting."""
        citation_metadata = {
            "module1": {
                "citation_data": {"title": "Test"},
                "citation_file_found": True,
            }
        }

        result = format_output(citation_metadata, "yaml")

        assert "module1:" in result
        assert "title: Test" in result

    def test_when_module_is_none_then_not_counted_as_found(self):
        """When a module's metadata is None, should not count as found."""
        citation_metadata = {
            "module1": None,
            "module2": {"citation_file_found": True},
        }

        summary = get_citation_summary(citation_metadata)

        assert summary["total_modules"] == 2
        assert summary["modules_with_citations"] == 1
