"""Tests for the new strategy-pattern-based validation system."""

import pytest
from omnibenchmark.benchmark.metadata import (
    ValidationResult,
    ValidationIssue,
    ValidationSeverity,
    StrictValidationStrategy,
    WarnValidationStrategy,
    ValidationContext,
    ValidationException,
    create_validation_context,
    validate_citation_cff_content,
    validate_citation_cff_data,
    validate_citation_authors,
    validate_citation_license,
    validate_omnibenchmark_yaml_content,
    validate_license_file_content,
    validate_license_consistency,
    validate_file_structure,
    detect_license_from_content,
    validate_module_files,
)


class TestValidationIssue:
    """Test ValidationIssue class."""

    def test_validation_issue_init(self):
        """Test ValidationIssue initialization."""
        issue = ValidationIssue(
            issue_type="test_issue",
            severity=ValidationSeverity.ERROR,
            path="test/path.txt",
            module_id="test_module",
            message="Test message",
            extra_data="test",
        )

        assert issue.issue_type == "test_issue"
        assert issue.severity == ValidationSeverity.ERROR
        assert issue.path == "test/path.txt"
        assert issue.module_id == "test_module"
        assert issue.message == "Test message"
        assert issue.context["extra_data"] == "test"
        assert issue.is_error is True
        assert issue.is_warning is False
        assert issue.msg == "Test message"  # backward compatibility

    def test_validation_issue_equality(self):
        """Test ValidationIssue equality and hashing."""
        issue1 = ValidationIssue(
            issue_type="test_issue",
            severity=ValidationSeverity.ERROR,
            path="test/path.txt",
            module_id="test_module",
            message="Test message",
        )
        issue2 = ValidationIssue(
            issue_type="test_issue",
            severity=ValidationSeverity.ERROR,
            path="test/path.txt",
            module_id="test_module",
            message="Test message",
        )
        issue3 = ValidationIssue(
            issue_type="different_issue",
            severity=ValidationSeverity.ERROR,
            path="test/path.txt",
            module_id="test_module",
            message="Test message",
        )

        assert issue1 == issue2
        assert issue1 != issue3
        assert hash(issue1) == hash(issue2)
        assert hash(issue1) != hash(issue3)


class TestValidationStrategies:
    """Test validation strategies."""

    def test_strict_validation_strategy(self):
        """Test StrictValidationStrategy."""
        strategy = StrictValidationStrategy()

        issue = strategy.handle_issue(
            issue_type="test_issue",
            path="test.txt",
            module_id="test_module",
            message="Test error",
        )

        assert issue.severity == ValidationSeverity.ERROR
        assert strategy.should_fail_fast() is True

    def test_warn_validation_strategy(self):
        """Test WarnValidationStrategy."""
        strategy = WarnValidationStrategy()

        issue = strategy.handle_issue(
            issue_type="test_issue",
            path="test.txt",
            module_id="test_module",
            message="Test warning",
        )

        assert issue.severity == ValidationSeverity.WARNING
        assert strategy.should_fail_fast() is False


class TestValidationResult:
    """Test ValidationResult class."""

    def test_validation_result_init(self):
        """Test ValidationResult initialization."""
        result = ValidationResult()

        assert len(result.issues) == 0
        assert len(result.errors) == 0
        assert len(result.warnings) == 0
        assert result.is_valid() is True
        assert result.has_warnings() is False

    def test_validation_result_add_issue(self):
        """Test adding issues to ValidationResult."""
        result = ValidationResult()

        error_issue = ValidationIssue(
            issue_type="test_error",
            severity=ValidationSeverity.ERROR,
            path="test.txt",
            module_id="test_module",
            message="Test error",
        )
        warning_issue = ValidationIssue(
            issue_type="test_warning",
            severity=ValidationSeverity.WARNING,
            path="test.txt",
            module_id="test_module",
            message="Test warning",
        )

        result.add_issue(error_issue)
        result.add_issue(warning_issue)

        assert len(result.issues) == 2
        assert len(result.errors) == 1
        assert len(result.warnings) == 1
        assert result.is_valid() is False
        assert result.has_warnings() is True

    def test_validation_result_deduplication(self):
        """Test that ValidationResult deduplicates identical issues."""
        result = ValidationResult()

        issue1 = ValidationIssue(
            issue_type="test_issue",
            severity=ValidationSeverity.ERROR,
            path="test.txt",
            module_id="test_module",
            message="Test message",
        )
        issue2 = ValidationIssue(
            issue_type="test_issue",
            severity=ValidationSeverity.ERROR,
            path="test.txt",
            module_id="test_module",
            message="Test message",
        )

        result.add_issue(issue1)
        result.add_issue(issue2)

        assert len(result.issues) == 1  # Deduplication worked

    def test_validation_result_get_issues_by_module(self):
        """Test grouping issues by module."""
        result = ValidationResult()

        issue1 = ValidationIssue(
            issue_type="test_issue1",
            severity=ValidationSeverity.ERROR,
            path="test1.txt",
            module_id="module1",
            message="Test message 1",
        )
        issue2 = ValidationIssue(
            issue_type="test_issue2",
            severity=ValidationSeverity.WARNING,
            path="test2.txt",
            module_id="module2",
            message="Test message 2",
        )
        issue3 = ValidationIssue(
            issue_type="test_issue3",
            severity=ValidationSeverity.ERROR,
            path="test3.txt",
            module_id="module1",
            message="Test message 3",
        )

        result.add_issue(issue1)
        result.add_issue(issue2)
        result.add_issue(issue3)

        by_module = result.get_issues_by_module()

        assert len(by_module) == 2
        assert len(by_module["module1"]) == 2
        assert len(by_module["module2"]) == 1


class TestValidationContext:
    """Test ValidationContext class."""

    def test_validation_context_strict_mode(self):
        """Test ValidationContext in strict mode."""
        ctx = ValidationContext(
            strategy=StrictValidationStrategy(),
            module_id="test_module",
            base_path="/test",
        )

        with pytest.raises(ValidationException) as exc_info:
            ctx.add_issue(
                issue_type="test_error", path="test.txt", message="Test error message"
            )

        assert "Test error message" in str(exc_info.value)
        assert len(exc_info.value.issues) == 1
        assert exc_info.value.issues[0].is_error

    def test_validation_context_warn_mode(self):
        """Test ValidationContext in warn mode."""
        ctx = ValidationContext(
            strategy=WarnValidationStrategy(),
            module_id="test_module",
            base_path="/test",
        )

        # Should not raise exception
        ctx.add_issue(issue_type="test_issue", path="test.txt", message="Test message")

        assert len(ctx.result.issues) == 1
        assert ctx.result.issues[0].is_warning

    def test_create_validation_context(self):
        """Test create_validation_context factory function."""
        strict_ctx = create_validation_context("module1", warn_mode=False)
        warn_ctx = create_validation_context("module2", warn_mode=True)

        assert isinstance(strict_ctx.strategy, StrictValidationStrategy)
        assert isinstance(warn_ctx.strategy, WarnValidationStrategy)
        assert strict_ctx.module_id == "module1"
        assert warn_ctx.module_id == "module2"


class TestCitationValidation:
    """Test citation validation functions."""

    def test_validate_citation_cff_content_valid(self):
        """Test validation with valid CITATION.cff content."""
        content = """
cff-version: 1.2.0
message: If you use this software, please cite it as below.
title: Test Software Package
authors:
  - family-names: Doe
    given-names: John
    orcid: https://orcid.org/0000-0000-0000-0000
  - family-names: Smith
    given-names: Jane
license: MIT
version: 1.0.0
date-released: 2023-01-01
"""
        ctx = create_validation_context("test_module", warn_mode=True)
        result, data = validate_citation_cff_content(content, ctx)

        assert result.is_valid()
        assert data is not None
        assert data["title"] == "Test Software Package"
        assert len(data["authors"]) == 2

    def test_validate_citation_cff_content_invalid_yaml(self):
        """Test validation with invalid YAML in CITATION.cff."""
        content = """
cff-version: 1.2.0
message: If you use this software, please cite it as below.
title: Test Software Package
invalid: yaml: structure: here
license: MIT
"""
        # Test strict mode - should raise exception (fail fast)
        ctx = create_validation_context("test_module", warn_mode=False)
        with pytest.raises(ValidationException) as exc_info:
            validate_citation_cff_content(content, ctx)

        assert len(exc_info.value.issues) == 1
        assert exc_info.value.issues[0].issue_type == "citation_invalid_yaml"

        # Test warn mode - should have warnings instead
        ctx_warn = create_validation_context("test_module", warn_mode=True)
        result_warn, data_warn = validate_citation_cff_content(content, ctx_warn)

        assert result_warn.is_valid()  # Valid because errors converted to warnings
        assert data_warn is None
        # Check that we have the issue in ctx_warn.result since result_warn might be empty
        assert len(ctx_warn.result.warnings) > 0
        assert any(
            issue.issue_type == "citation_invalid_yaml"
            for issue in ctx_warn.result.warnings
        )

    def test_validate_citation_cff_data_missing_required_fields(self):
        """Test validation with missing required fields."""
        data = {
            "cff-version": "1.2.0",
            # Missing: message, title, authors
            "license": "MIT",
        }
        # Test strict mode - should raise exception on first missing field
        ctx = create_validation_context("test_module", warn_mode=False)
        with pytest.raises(ValidationException) as exc_info:
            validate_citation_cff_data(data, ctx)

        assert len(exc_info.value.issues) == 1
        assert exc_info.value.issues[0].issue_type == "citation_missing_required_field"

    def test_validate_citation_authors_valid(self):
        """Test validation with valid authors."""
        authors = [
            {"family-names": "Doe", "given-names": "John"},
            {"family-names": "Smith", "given-names": "Jane"},
        ]
        ctx = create_validation_context("test_module", warn_mode=True)
        result = validate_citation_authors(authors, ctx)

        assert result.is_valid()

    def test_validate_citation_authors_missing_family_names(self):
        """Test validation with authors missing family names."""
        authors = [
            {"given-names": "John"},  # Missing family-names
            {"family-names": "Smith", "given-names": "Jane"},
        ]
        # Test strict mode - should raise exception on first error
        ctx = create_validation_context("test_module", warn_mode=False)
        with pytest.raises(ValidationException) as exc_info:
            validate_citation_authors(authors, ctx)

        assert len(exc_info.value.issues) == 1
        assert (
            exc_info.value.issues[0].issue_type == "citation_author_missing_family_name"
        )

    def test_validate_citation_authors_missing_given_names(self):
        """Test validation with authors missing given names (warning only)."""
        authors = [
            {"family-names": "Doe"},  # Missing given-names
            {"family-names": "Smith", "given-names": "Jane"},
        ]
        ctx = create_validation_context(
            "test_module", warn_mode=False
        )  # Even in strict mode
        result = validate_citation_authors(authors, ctx)

        assert result.is_valid()  # Should be valid since given-names is optional
        assert result.has_warnings()  # But should have warnings
        assert any(
            issue.issue_type == "citation_author_missing_given_name"
            for issue in result.warnings
        )

    def test_validate_citation_license_valid(self):
        """Test validation with valid license."""
        ctx = create_validation_context("test_module", warn_mode=True)
        result = validate_citation_license("MIT", ctx)

        assert result.is_valid()

    def test_validate_citation_license_invalid_spdx(self):
        """Test validation with invalid SPDX license."""
        # Test strict mode - should raise exception
        ctx = create_validation_context("test_module", warn_mode=False)
        with pytest.raises(ValidationException) as exc_info:
            validate_citation_license("INVALID-LICENSE", ctx)

        assert len(exc_info.value.issues) == 1
        assert exc_info.value.issues[0].issue_type == "citation_invalid_license"


class TestOmnibenchmarkValidation:
    """Test omnibenchmark.yaml validation functions."""

    def test_validate_omnibenchmark_yaml_content_valid(self):
        """Test validation with valid omnibenchmark.yaml content."""
        content = """
name: test-software-package
version: 1.0.0
description: A comprehensive test software package
type: method
stage: preprocessing
"""
        ctx = create_validation_context("test_module", warn_mode=True)
        result, data = validate_omnibenchmark_yaml_content(content, ctx)

        assert result.is_valid()
        assert data is not None
        assert data["name"] == "test-software-package"

    def test_validate_omnibenchmark_yaml_content_invalid_yaml(self):
        """Test validation with invalid YAML in omnibenchmark.yaml."""
        content = """
name: test-software-package
version: 1.0.0
invalid: yaml: structure: here
type: method
"""
        ctx = create_validation_context("test_module", warn_mode=True)
        result, data = validate_omnibenchmark_yaml_content(content, ctx)

        assert data is None
        assert result.has_warnings()  # Should be warnings, not errors
        assert any(
            issue.issue_type == "omnibenchmark_yaml_invalid_yaml"
            for issue in result.warnings
        )


class TestLicenseValidation:
    """Test license validation functions."""

    def test_validate_license_file_content_valid(self):
        """Test validation with valid LICENSE content."""
        content = """
MIT License

Copyright (c) 2023 John Doe

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
        ctx = create_validation_context("test_module", warn_mode=True)
        result = validate_license_file_content(content, ctx)

        assert result.is_valid()

    def test_validate_license_file_content_empty(self):
        """Test validation with empty LICENSE file."""
        # Test strict mode - should raise exception
        ctx = create_validation_context("test_module", warn_mode=False)
        with pytest.raises(ValidationException) as exc_info:
            validate_license_file_content("", ctx)

        assert len(exc_info.value.issues) == 1
        assert exc_info.value.issues[0].issue_type == "license_file_empty"

    def test_validate_license_consistency_both_present(self):
        """Test license consistency when both citation and file are present."""
        citation_license = "MIT"
        license_content = "MIT License\n\nPermission is hereby granted..."

        ctx = create_validation_context("test_module", warn_mode=True)
        result = validate_license_consistency(citation_license, license_content, ctx)

        assert result.is_valid()

    def test_validate_license_consistency_both_missing(self):
        """Test license consistency when both are missing."""
        # Test strict mode - should raise exception
        ctx = create_validation_context("test_module", warn_mode=False)
        with pytest.raises(ValidationException) as exc_info:
            validate_license_consistency(None, None, ctx)

        assert len(exc_info.value.issues) == 1
        assert exc_info.value.issues[0].issue_type == "no_license_information"

    def test_detect_license_from_content_mit(self):
        """Test license detection for MIT license."""
        content = "MIT License\n\nPermission is hereby granted..."
        license_type = detect_license_from_content(content)
        assert license_type == "MIT"

    def test_detect_license_from_content_apache(self):
        """Test license detection for Apache license."""
        content = "Apache License Version 2.0\n\nLicensed under the Apache License..."
        license_type = detect_license_from_content(content)
        assert license_type == "Apache-2.0"

    def test_detect_license_from_content_gpl3(self):
        """Test license detection for GPL v3 license."""
        content = (
            "GNU GENERAL PUBLIC LICENSE Version 3\n\nThis program is free software..."
        )
        license_type = detect_license_from_content(content)
        assert license_type == "GPL-3.0"

    def test_detect_license_from_content_unknown(self):
        """Test license detection for unknown license."""
        content = "Some custom license text that doesn't match known patterns"
        license_type = detect_license_from_content(content)
        assert license_type is None


class TestFileStructureValidation:
    """Test file structure validation."""

    def test_validate_file_structure_all_present(self):
        """Test validation when all expected files are present."""
        files_present = {
            "CITATION.cff": True,
            "LICENSE": True,
            "omnibenchmark.yaml": True,
        }
        ctx = create_validation_context("test_module", warn_mode=True)
        result = validate_file_structure(files_present, ctx)

        assert result.is_valid()

    def test_validate_file_structure_missing_citation(self):
        """Test validation when CITATION.cff is missing."""
        files_present = {
            "CITATION.cff": False,
            "LICENSE": True,
            "omnibenchmark.yaml": True,
        }
        # Test strict mode - should raise exception
        ctx = create_validation_context("test_module", warn_mode=False)
        with pytest.raises(ValidationException) as exc_info:
            validate_file_structure(files_present, ctx)

        assert len(exc_info.value.issues) == 1
        assert exc_info.value.issues[0].issue_type == "citation_missing"

    def test_validate_file_structure_missing_license(self):
        """Test validation when LICENSE file is missing."""
        files_present = {
            "CITATION.cff": True,
            "LICENSE": False,
            "LICENSE.txt": False,
            "LICENSE.md": False,
            "omnibenchmark.yaml": True,
        }
        ctx = create_validation_context("test_module", warn_mode=True)
        result = validate_file_structure(files_present, ctx)

        assert result.is_valid()  # Should be valid
        assert result.has_warnings()  # But have warnings
        assert any(issue.issue_type == "no_license_file" for issue in result.warnings)


class TestIntegrationValidation:
    """Integration tests for comprehensive validation."""

    def test_validate_module_files_all_valid(self):
        """Test comprehensive validation with all valid files."""
        citation_content = """
cff-version: 1.2.0
message: If you use this software, please cite it as below.
title: Test Software Package
authors:
  - family-names: Doe
    given-names: John
license: MIT
"""

        license_content = """
MIT License

Copyright (c) 2023 John Doe

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

        omnibenchmark_content = """
name: test-software-package
version: 1.0.0
description: A comprehensive test software package
type: method
stage: preprocessing
"""

        files_present = {
            "CITATION.cff": True,
            "LICENSE": True,
            "omnibenchmark.yaml": True,
        }

        result = validate_module_files(
            module_id="test_module",
            citation_content=citation_content,
            license_content=license_content,
            omnibenchmark_content=omnibenchmark_content,
            files_present=files_present,
            warn_mode=True,
        )

        assert result.is_valid()

    def test_validate_module_files_with_errors_and_warnings(self):
        """Test comprehensive validation with mixed errors and warnings."""
        citation_content = """
cff-version: 1.2.0
# Missing required fields: message, title, authors
license: INVALID-LICENSE
"""

        files_present = {
            "CITATION.cff": True,
            "LICENSE": False,  # Missing LICENSE file
            "omnibenchmark.yaml": False,  # Missing omnibenchmark.yaml
        }

        # Test in strict mode first - should have errors
        result_strict = validate_module_files(
            module_id="test_module",
            citation_content=citation_content,
            files_present=files_present,
            warn_mode=False,
        )

        assert not result_strict.is_valid()
        assert len(result_strict.errors) > 0

        # Test in warn mode - should convert most errors to warnings
        # BUT omnibenchmark.yaml missing is ALWAYS an error
        result_warn = validate_module_files(
            module_id="test_module",
            citation_content=citation_content,
            files_present=files_present,
            warn_mode=True,
        )

        assert (
            not result_warn.is_valid()
        )  # Invalid because omnibenchmark.yaml is always an error
        assert result_warn.has_warnings()
        # Should have at least one error (omnibenchmark.yaml missing)
        assert len(result_warn.errors) >= 1
        assert any(
            e.issue_type == "omnibenchmark_yaml_missing" for e in result_warn.errors
        )

        # Should have issues for missing required fields, invalid license, missing files
        issue_types = {issue.issue_type for issue in result_warn.issues}
        expected_types = {
            "citation_missing_required_field",
            "citation_invalid_license",
            "no_license_file",
            "omnibenchmark_yaml_missing",
        }

        # Check that we have some of the expected issue types
        assert len(issue_types.intersection(expected_types)) > 0

    def test_validate_module_files_strict_mode_fails_fast(self):
        """Test that strict mode fails fast on first error."""
        citation_content = """
cff-version: 1.2.0
# Missing required fields will cause errors
"""

        files_present = {
            "CITATION.cff": True,
            "LICENSE": True,
            "omnibenchmark.yaml": True,
        }

        # In strict mode, should fail fast
        result = validate_module_files(
            module_id="test_module",
            citation_content=citation_content,
            files_present=files_present,
            warn_mode=False,  # Strict mode
        )

        # Should have collected some errors even if it failed fast
        assert not result.is_valid()
        assert len(result.errors) > 0
