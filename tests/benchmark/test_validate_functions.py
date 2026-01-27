"""Unit tests for benchmark validation functions."""

import pytest

from omnibenchmark.benchmark.validate import (
    ValidationResult,
    _detect_license_from_content,
    _validate_license_consistency,
    _validate_citation_cff,
    _validate_license_file,
    _validate_omnibenchmark_yaml,
    get_validation_summary,
    format_validation_results,
)


class TestValidationResult:
    """Unit tests for ValidationResult class."""

    @pytest.mark.short
    def test_init_creates_empty_result(self):
        """Test that ValidationResult initializes with empty lists."""
        result = ValidationResult("test_module")

        assert result.module_id == "test_module"
        assert result.validation_errors == []
        assert result.validation_warnings == []
        assert result.repository_url is None
        assert result.commit_hash is None

    @pytest.mark.short
    def test_add_error(self):
        """Test adding validation errors."""
        result = ValidationResult("test_module")

        result.add_error("Test error message")

        assert len(result.validation_errors) == 1
        assert "Test error message" in result.validation_errors

    @pytest.mark.short
    def test_add_warning(self):
        """Test adding validation warnings."""
        result = ValidationResult("test_module")

        result.add_warning("Test warning message")

        assert len(result.validation_warnings) == 1
        assert "Test warning message" in result.validation_warnings

    @pytest.mark.short
    def test_is_valid_with_no_errors(self):
        """Test that result is valid when no errors and all checks pass."""
        result = ValidationResult("test_module")
        result.repository_url = "https://github.com/test/repo"
        result.commit_hash = "abc123"
        result.local_repo_exists = True
        result.citation_file_exists = True
        result.citation_file_valid = True
        result.citation_has_license = True
        result.citation_has_authors = True

        assert result.is_valid()

    @pytest.mark.short
    def test_is_valid_false_with_errors(self):
        """Test that result is invalid when errors exist."""
        result = ValidationResult("test_module")
        result.repository_url = "https://github.com/test/repo"
        result.commit_hash = "abc123"
        result.local_repo_exists = True
        result.citation_file_exists = True
        result.citation_file_valid = True
        result.citation_has_license = True
        result.citation_has_authors = True
        result.add_error("Test error")

        assert not result.is_valid()


class TestDetectLicenseFromContent:
    """Unit tests for license detection."""

    @pytest.mark.short
    def test_detect_mit_license(self):
        """Test detection of MIT license."""
        content = """MIT License

Copyright (c) 2023 Test Author

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction."""

        detected = _detect_license_from_content(content)

        assert detected == "MIT"

    @pytest.mark.short
    def test_detect_apache_license(self):
        """Test detection of Apache 2.0 license."""
        content = """Apache License
Version 2.0, January 2004
http://www.apache.org/licenses/"""

        detected = _detect_license_from_content(content)

        assert detected == "Apache-2.0"

    @pytest.mark.short
    def test_detect_gpl3_license(self):
        """Test detection of GPL-3.0 license."""
        content = """GNU GENERAL PUBLIC LICENSE
Version 3, 29 June 2007"""

        detected = _detect_license_from_content(content)

        assert detected == "GPL-3.0"

    @pytest.mark.short
    def test_detect_unknown_license(self):
        """Test that unknown licenses return None."""
        content = "Some random text without license keywords"

        detected = _detect_license_from_content(content)

        assert detected is None


class TestValidateLicenseConsistency:
    """Unit tests for license consistency validation."""

    @pytest.mark.short
    def test_consistent_licenses_no_warning(self):
        """Test that consistent licenses don't generate warnings."""
        result = ValidationResult("test_module")
        result.citation_license = "MIT"
        result.citation_has_license = True
        result.license_file_exists = True
        result.license_file_content = """MIT License

Copyright (c) 2023 Test Author

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software."""

        _validate_license_consistency(result)

        assert len(result.validation_warnings) == 0
        assert len(result.validation_errors) == 0
        assert result.license_matches_citation

    @pytest.mark.short
    def test_inconsistent_licenses_error(self):
        """Test that inconsistent licenses generate errors."""
        result = ValidationResult("test_module")
        result.citation_license = "MIT"
        result.citation_has_license = True
        result.license_file_exists = True
        result.license_file_content = "Apache License\nVersion 2.0"

        _validate_license_consistency(result)

        assert len(result.validation_errors) > 0
        assert not result.license_matches_citation
        assert any("mismatch" in e.lower() for e in result.validation_errors)

    @pytest.mark.short
    def test_undetectable_license_warning(self):
        """Test that undetectable license generates warning."""
        result = ValidationResult("test_module")
        result.citation_license = "MIT"
        result.citation_has_license = True
        result.license_file_exists = True
        result.license_file_content = "Some custom license text"

        _validate_license_consistency(result)

        assert len(result.validation_warnings) > 0
        assert any("could not" in w.lower() for w in result.validation_warnings)

    @pytest.mark.short
    def test_no_citation_license_skips_check(self):
        """Test that missing citation license skips consistency check."""
        result = ValidationResult("test_module")
        result.citation_has_license = False
        result.license_file_exists = True
        result.license_file_content = "MIT License"

        _validate_license_consistency(result)

        # No consistency check if citation has no license
        assert len(result.validation_errors) == 0
        assert len(result.validation_warnings) == 0


class TestValidateCitationCff:
    """Unit tests for CITATION.cff validation."""

    @pytest.mark.short
    def test_validate_citation_with_valid_file(self, tmp_path):
        """Test validation of valid CITATION.cff file."""
        citation_file = tmp_path / "CITATION.cff"
        citation_file.write_text("""cff-version: 1.2.0
message: Please cite this software
title: Test Software
authors:
  - family-names: Doe
    given-names: John
license: MIT
""")

        result = ValidationResult("test_module")
        _validate_citation_cff(tmp_path, result)

        assert result.citation_file_exists
        assert result.citation_file_valid
        assert result.citation_has_authors
        assert result.citation_has_license
        assert result.citation_license == "MIT"

    @pytest.mark.short
    def test_validate_citation_missing_file(self, tmp_path):
        """Test validation when CITATION.cff is missing."""
        result = ValidationResult("test_module")
        _validate_citation_cff(tmp_path, result)

        assert not result.citation_file_exists
        assert len(result.validation_errors) > 0


class TestValidateLicenseFile:
    """Unit tests for LICENSE file validation."""

    @pytest.mark.short
    def test_validate_license_file_exists(self, tmp_path):
        """Test validation when LICENSE file exists."""
        license_file = tmp_path / "LICENSE"
        license_file.write_text("MIT License\n\nCopyright (c) 2023")

        result = ValidationResult("test_module")
        _validate_license_file(tmp_path, result)

        assert result.license_file_exists
        assert result.license_file_content is not None

    @pytest.mark.short
    def test_validate_license_file_missing_with_citation_license(self, tmp_path):
        """Test validation when LICENSE file is missing but citation has license."""
        result = ValidationResult("test_module")
        result.citation_has_license = True  # Citation has license
        _validate_license_file(tmp_path, result)

        assert not result.license_file_exists
        assert len(result.validation_warnings) > 0

    @pytest.mark.short
    def test_validate_license_file_missing_no_citation_license(self, tmp_path):
        """Test validation when LICENSE file and citation license are both missing."""
        result = ValidationResult("test_module")
        result.citation_has_license = False  # No license in citation
        _validate_license_file(tmp_path, result)

        assert not result.license_file_exists
        assert len(result.validation_errors) > 0


class TestValidateOmnibenchmarkYaml:
    """Unit tests for omnibenchmark.yaml validation."""

    @pytest.mark.short
    def test_validate_omnibenchmark_yaml_exists(self, tmp_path):
        """Test validation when omnibenchmark.yaml exists."""
        yaml_file = tmp_path / "omnibenchmark.yaml"
        yaml_file.write_text("name: test_module\nversion: '1.0.0'\n")

        result = ValidationResult("test_module")
        _validate_omnibenchmark_yaml(tmp_path, result)

        assert result.omnibenchmark_yaml_exists
        assert result.omnibenchmark_yaml_valid

    @pytest.mark.short
    def test_validate_omnibenchmark_yaml_missing(self, tmp_path):
        """Test validation when omnibenchmark.yaml is missing."""
        result = ValidationResult("test_module")
        _validate_omnibenchmark_yaml(tmp_path, result)

        assert not result.omnibenchmark_yaml_exists


class TestGetValidationSummary:
    """Unit tests for validation summary generation."""

    @pytest.mark.short
    def test_summary_all_valid(self):
        """Test summary when all modules are valid."""
        results = {
            "module1": ValidationResult("module1"),
            "module2": ValidationResult("module2"),
        }
        # Make them valid
        for result in results.values():
            result.repository_url = "https://github.com/test/repo"
            result.commit_hash = "abc123"
            result.local_repo_exists = True
            result.citation_file_exists = True
            result.citation_file_valid = True
            result.citation_has_license = True
            result.citation_has_authors = True

        summary = get_validation_summary(results)

        assert summary["total_modules"] == 2
        assert summary["valid_modules"] == 2
        assert summary["validation_success_rate"] == 1.0

    @pytest.mark.short
    def test_summary_with_errors(self):
        """Test summary when some modules have errors."""
        results = {
            "module1": ValidationResult("module1"),
            "module2": ValidationResult("module2"),
        }
        # Make module2 valid
        results["module2"].repository_url = "https://github.com/test/repo"
        results["module2"].commit_hash = "abc123"
        results["module2"].local_repo_exists = True
        results["module2"].citation_file_exists = True
        results["module2"].citation_file_valid = True
        results["module2"].citation_has_license = True
        results["module2"].citation_has_authors = True

        # module1 has error so is invalid
        results["module1"].add_error("Test error")

        summary = get_validation_summary(results)

        assert summary["total_modules"] == 2
        assert summary["valid_modules"] == 1


class TestFormatValidationResults:
    """Unit tests for formatting validation results."""

    @pytest.mark.short
    def test_format_results_json(self):
        """Test JSON formatting of validation results."""
        results = {
            "module1": ValidationResult("module1"),
        }
        results["module1"].add_error("Test error")

        output = format_validation_results(results, format_type="json")

        assert '"module1"' in output
        assert "Test error" in output

    @pytest.mark.short
    def test_format_results_text(self):
        """Test text formatting of validation results."""
        results = {
            "module1": ValidationResult("module1"),
        }
        results["module1"].add_error("Test error")

        output = format_validation_results(results, format_type="text")

        assert "module1" in output
        assert "Total modules" in output
