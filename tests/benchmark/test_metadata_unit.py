"""Unit tests for metadata validation functions with minimal mocking."""

from omnibenchmark.benchmark.metadata import (
    validate_citation_cff_content,
    validate_license_file_content,
    validate_omnibenchmark_yaml_content,
    validate_file_structure,
    validate_license_consistency,
    validate_module_files,
    ValidationResult,
    ValidationIssue,
    ValidationSeverity,
    ValidationException,
    create_validation_context,
)


class TestValidateCitationCffContent:
    """Unit tests for CITATION.cff content validation."""

    def test_valid_citation_cff(self):
        """Test validation of valid CITATION.cff content."""
        citation_content = """
cff-version: 1.2.0
message: If you use this software, please cite it as below.
title: Test Module
authors:
  - family-names: Doe
    given-names: John
    email: john.doe@example.com
license: MIT
version: 1.0.0
date-released: 2023-01-01
"""
        ctx = create_validation_context("test_module", warn_mode=False)
        result, data = validate_citation_cff_content(
            citation_content, ctx, "CITATION.cff"
        )

        assert result.is_valid()
        assert len(result.errors) == 0
        assert isinstance(data, dict)
        assert data["title"] == "Test Module"
        assert data["license"] == "MIT"

    def test_missing_required_fields(self):
        """Test validation fails when required fields are missing."""
        citation_content = """
cff-version: 1.2.0
title: Test Module
"""
        ctx = create_validation_context("test_module", warn_mode=True)
        result, data = validate_citation_cff_content(
            citation_content, ctx, "CITATION.cff"
        )

        # Check the context result, not the function return result
        assert ctx.result.is_valid()  # Valid in warn mode
        assert ctx.result.has_warnings()

        warning_messages = [issue.message for issue in ctx.result.warnings]
        assert any("message" in msg.lower() for msg in warning_messages)
        assert any("author" in msg.lower() for msg in warning_messages)

    def test_invalid_yaml_syntax(self):
        """Test validation handles invalid YAML syntax."""
        citation_content = """
cff-version: 1.2.0
title: Test Module
invalid: yaml: [
"""
        ctx = create_validation_context("test_module", warn_mode=True)
        result, data = validate_citation_cff_content(
            citation_content, ctx, "CITATION.cff"
        )

        assert ctx.result.is_valid()  # Valid in warn mode
        assert ctx.result.has_warnings()
        assert data is None

    def test_warn_mode_converts_errors_to_warnings(self):
        """Test that warn mode converts errors to warnings."""
        citation_content = """
cff-version: 1.2.0
title: Test Module
"""
        ctx = create_validation_context("test_module", warn_mode=True)
        result, data = validate_citation_cff_content(
            citation_content, ctx, "CITATION.cff"
        )

        assert result.is_valid()  # Valid in warn mode
        assert result.has_warnings()
        assert len(result.errors) == 0


class TestValidateLicenseFileContent:
    """Unit tests for LICENSE file content validation."""

    def test_valid_license_file(self):
        """Test validation of valid LICENSE file."""
        license_content = """MIT License

Copyright (c) 2023 Test Author

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
"""
        ctx = create_validation_context("test_module", warn_mode=False)
        result = validate_license_file_content(license_content, ctx, "LICENSE")

        assert result.is_valid()
        assert len(result.errors) == 0

    def test_empty_license_file(self):
        """Test validation fails for empty LICENSE file."""
        ctx = create_validation_context("test_module", warn_mode=True)
        validate_license_file_content("", ctx, "LICENSE")

        assert ctx.result.is_valid()  # Valid in warn mode
        assert ctx.result.has_warnings()

    def test_very_short_license_file(self):
        """Test validation fails for very short LICENSE file."""
        ctx = create_validation_context("test_module", warn_mode=True)
        validate_license_file_content("MIT", ctx, "LICENSE")

        assert ctx.result.is_valid()  # Valid in warn mode
        assert ctx.result.has_warnings()


class TestValidateOmnibenchmarkYamlContent:
    """Unit tests for omnibenchmark.yaml content validation."""

    def test_valid_omnibenchmark_yaml(self):
        """Test validation of valid omnibenchmark.yaml."""
        yaml_content = """
name: test_module
version: "1.0.0"
description: A test module for benchmarking
"""
        ctx = create_validation_context("test_module", warn_mode=False)
        result, data = validate_omnibenchmark_yaml_content(
            yaml_content, ctx, "omnibenchmark.yaml"
        )

        assert result.is_valid()
        assert len(result.errors) == 0
        assert isinstance(data, dict)
        assert data["name"] == "test_module"
        assert data["version"] == "1.0.0"

    def test_missing_name_field(self):
        """Test validation fails when name field is missing."""
        yaml_content = """
version: "1.0.0"
description: A test module
"""
        ctx = create_validation_context("test_module", warn_mode=True)
        result, data = validate_omnibenchmark_yaml_content(
            yaml_content, ctx, "omnibenchmark.yaml"
        )

        assert ctx.result.is_valid()  # Valid in warn mode
        assert ctx.result.has_warnings()

    def test_invalid_version_format(self):
        """Test validation fails for invalid version format."""
        yaml_content = """
name: test_module
version: 1.0  # Should be string
"""
        ctx = create_validation_context("test_module", warn_mode=True)
        result, data = validate_omnibenchmark_yaml_content(
            yaml_content, ctx, "omnibenchmark.yaml"
        )

        assert ctx.result.is_valid()  # Valid in warn mode
        assert ctx.result.has_warnings()

    def test_empty_omnibenchmark_yaml(self):
        """Test validation handles empty omnibenchmark.yaml."""
        ctx = create_validation_context("test_module", warn_mode=True)
        result, data = validate_omnibenchmark_yaml_content(
            "", ctx, "omnibenchmark.yaml"
        )

        assert ctx.result.is_valid()  # Valid in warn mode
        assert ctx.result.has_warnings()
        assert data is None


class TestValidateFileStructure:
    """Unit tests for file structure validation."""

    def test_all_files_present(self):
        """Test validation passes when all files are present."""
        files_present = {
            "CITATION.cff": True,
            "LICENSE": True,
            "omnibenchmark.yaml": True,
        }
        ctx = create_validation_context("test_module", warn_mode=False)
        result = validate_file_structure(files_present, ctx, citation_has_license=False)

        assert result.is_valid()
        assert not result.has_warnings()

    def test_missing_citation_file(self):
        """Test validation fails when CITATION.cff is missing."""
        files_present = {
            "CITATION.cff": False,
            "LICENSE": True,
            "omnibenchmark.yaml": True,
        }
        ctx = create_validation_context("test_module", warn_mode=True)
        validate_file_structure(files_present, ctx, citation_has_license=False)

        assert ctx.result.is_valid()  # Valid in warn mode
        assert ctx.result.has_warnings()

    def test_missing_license_file_with_citation_license(self):
        """Test no warning for missing LICENSE when citation has license."""
        files_present = {
            "CITATION.cff": True,
            "LICENSE": False,
            "omnibenchmark.yaml": True,
        }
        ctx = create_validation_context("test_module", warn_mode=False)
        result = validate_file_structure(files_present, ctx, citation_has_license=True)

        assert result.is_valid()
        # Should not warn about missing LICENSE since citation has license

    def test_missing_license_file_without_citation_license(self):
        """Test warning for missing LICENSE when citation has no license."""
        files_present = {
            "CITATION.cff": True,
            "LICENSE": False,
            "omnibenchmark.yaml": True,
        }
        ctx = create_validation_context("test_module", warn_mode=False)
        result = validate_file_structure(files_present, ctx, citation_has_license=False)

        assert result.is_valid()
        assert result.has_warnings()

    def test_missing_omnibenchmark_yaml_with_legacy_config(self):
        """Test warning for legacy config.cfg instead of omnibenchmark.yaml."""
        files_present = {
            "CITATION.cff": True,
            "LICENSE": True,
            "omnibenchmark.yaml": False,
            "config.cfg": True,
        }
        ctx = create_validation_context("test_module", warn_mode=False)
        result = validate_file_structure(files_present, ctx, citation_has_license=False)

        assert result.is_valid()
        assert result.has_warnings()

        warning_messages = [issue.message for issue in result.warnings]
        assert any("legacy" in msg.lower() for msg in warning_messages)


class TestValidateLicenseConsistency:
    """Unit tests for license consistency validation."""

    def test_consistent_licenses(self):
        """Test validation passes when licenses are consistent."""
        license_content = """MIT License

Copyright (c) 2023 Test Author
"""
        ctx = create_validation_context("test_module", warn_mode=False)
        result = validate_license_consistency("MIT", license_content, ctx)

        assert result.is_valid()
        assert not result.has_warnings()

    def test_inconsistent_licenses(self):
        """Test validation warns when licenses are inconsistent."""
        license_content = """Apache License
Version 2.0, January 2004
"""
        ctx = create_validation_context("test_module", warn_mode=False)
        result = validate_license_consistency("MIT", license_content, ctx)

        assert result.is_valid()
        assert result.has_warnings()

    def test_no_citation_license(self):
        """Test validation warns when citation has no license but LICENSE file exists."""
        license_content = """MIT License

Copyright (c) 2023 Test Author
"""
        ctx = create_validation_context("test_module", warn_mode=False)
        result = validate_license_consistency(None, license_content, ctx)

        assert result.is_valid()
        assert result.has_warnings()  # Should warn about missing license in citation


class TestValidateModuleFiles:
    """Integration tests for complete module validation."""

    def test_complete_valid_module(self):
        """Test validation of completely valid module."""
        citation_content = """
cff-version: 1.2.0
message: If you use this software, please cite it as below.
title: Test Module
authors:
  - family-names: Doe
    given-names: John
license: MIT
"""
        license_content = """MIT License

Copyright (c) 2023 Test Author

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software").
"""
        omnibenchmark_content = """
name: test_module
version: "1.0.0"
entrypoints:
  default: script.py
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
            warn_mode=False,
        )

        assert result.is_valid()
        assert len(result.errors) == 0
        assert not result.has_warnings()

    def test_module_with_errors(self):
        """Test validation of module with errors."""
        citation_content = """
cff-version: 1.2.0
title: Test Module
# Missing required fields: message, authors
"""
        files_present = {
            "CITATION.cff": True,
            "LICENSE": False,
            "omnibenchmark.yaml": False,
        }

        result = validate_module_files(
            module_id="test_module",
            citation_content=citation_content,
            license_content=None,
            omnibenchmark_content=None,
            files_present=files_present,
            warn_mode=True,
        )

        assert (
            not result.is_valid()
        )  # Invalid because omnibenchmark.yaml is always an error
        assert result.has_warnings()
        # Should have error for missing omnibenchmark.yaml
        assert any(e.issue_type == "omnibenchmark_yaml_missing" for e in result.errors)

    def test_warn_mode_behavior(self):
        """Test validation in warn mode."""
        citation_content = """
cff-version: 1.2.0
title: Test Module
# Missing required fields but in warn mode
"""
        files_present = {
            "CITATION.cff": True,
            "LICENSE": True,
            "omnibenchmark.yaml": True,
        }

        result = validate_module_files(
            module_id="test_module",
            citation_content=citation_content,
            license_content="MIT License",
            omnibenchmark_content="name: test\nversion: '1.0'",
            files_present=files_present,
            warn_mode=True,
        )

        assert result.is_valid()  # Valid in warn mode
        assert result.has_warnings()


class TestValidationContext:
    """Unit tests for ValidationContext."""

    def test_create_validation_context(self):
        """Test creation of validation context."""
        ctx = create_validation_context(
            "test_module", warn_mode=True, base_path="/test"
        )

        assert ctx.module_id == "test_module"
        assert ctx.warn_mode is True
        assert ctx.base_path == "/test"
        assert isinstance(ctx.result, ValidationResult)

    def test_add_issue(self):
        """Test adding issues to context."""
        ctx = create_validation_context("test_module", warn_mode=True)

        ctx.add_issue(
            issue_type="test_error",
            path="test.txt",
            message="Test error message",
        )

        # In warn mode, errors are converted to warnings
        assert len(ctx.result.warnings) > 0
        assert len(ctx.result.warnings) == 1
        assert ctx.result.warnings[0].message == "Test error message"

    def test_merge_result(self):
        """Test merging validation results."""
        ctx = create_validation_context("test_module", warn_mode=True)

        other_result = ValidationResult()
        other_result.add_error("Another error")

        ctx.merge_result(other_result)

        assert len(ctx.result.errors) > 0
        assert len(ctx.result.errors) == 1


class TestValidationResult:
    """Unit tests for ValidationResult."""

    def test_validation_result_creation(self):
        """Test creation and basic properties of ValidationResult."""
        result = ValidationResult()

        assert result.is_valid()
        assert len(result.errors) == 0
        assert not result.has_warnings()
        assert len(result.errors) == 0
        assert len(result.warnings) == 0

    def test_add_error(self):
        """Test adding errors to validation result."""
        result = ValidationResult()

        # Create and add an error issue
        from omnibenchmark.benchmark.metadata import ValidationIssue, ValidationSeverity

        error_issue = ValidationIssue(
            issue_type="test_error",
            severity=ValidationSeverity.ERROR,
            path="test.txt",
            message="Test error",
            module_id="test_module",
        )
        result.add_issue(error_issue)

        assert not result.is_valid()
        assert len(result.errors) > 0
        assert len(result.errors) == 1
        assert result.errors[0].message == "Test error"

    def test_add_warning(self):
        """Test adding warnings to validation result."""
        result = ValidationResult()

        # Create and add a warning issue
        from omnibenchmark.benchmark.metadata import ValidationIssue, ValidationSeverity

        warning_issue = ValidationIssue(
            issue_type="test_warning",
            severity=ValidationSeverity.WARNING,
            path="test.txt",
            message="Test warning",
            module_id="test_module",
        )
        result.add_issue(warning_issue)

        assert result.is_valid()  # Warnings don't make it invalid
        assert result.has_warnings()
        assert len(result.warnings) == 1
        assert result.warnings[0].message == "Test warning"

    def test_add_issue(self):
        """Test adding ValidationIssue objects."""
        result = ValidationResult()

        error_issue = ValidationIssue(
            issue_type="test_error",
            severity=ValidationSeverity.ERROR,
            path="test.txt",
            message="Test error",
            module_id="test_module",
        )

        warning_issue = ValidationIssue(
            issue_type="test_warning",
            severity=ValidationSeverity.WARNING,
            path="test.txt",
            message="Test warning",
            module_id="test_module",
        )

        result.add_issue(error_issue)
        result.add_issue(warning_issue)

        assert not result.is_valid()
        assert len(result.errors) > 0
        assert result.has_warnings()
        assert len(result.errors) == 1
        assert len(result.warnings) == 1


class TestValidationException:
    """Unit tests for ValidationException."""

    def test_validation_exception_creation(self):
        """Test creation of ValidationException."""
        issue = ValidationIssue(
            issue_type="test_error",
            severity=ValidationSeverity.ERROR,
            path="test.txt",
            message="Test error",
            module_id="test_module",
        )

        exc = ValidationException("Test exception", [issue])

        assert str(exc) == "Test exception"
        assert len(exc.issues) == 1
        assert exc.issues[0] == issue

    def test_validation_exception_from_result(self):
        """Test creating ValidationException from ValidationResult."""
        result = ValidationResult()
        result.add_error("Test error")
        result.add_warning("Test warning")

        exc = ValidationException("Test", result.errors + result.warnings)

        assert len(exc.issues) == 2
        assert any(issue.severity == ValidationSeverity.ERROR for issue in exc.issues)
        assert any(issue.severity == ValidationSeverity.WARNING for issue in exc.issues)
