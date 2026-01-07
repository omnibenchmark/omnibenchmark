import pytest
from omnibenchmark.benchmark.metadata import (
    ValidationResult,
    ValidationIssue,
    ValidationSeverity,
    ValidationException,
    create_validation_context,
    validate_citation_cff_content,
    validate_license_consistency,
    validate_file_structure,
    validate_module_files,
)


@pytest.mark.short
def test_validate_complete_valid_module(tmp_path):
    """Test validation of a complete valid module using temporary files."""
    # Create valid CITATION.cff
    citation_content = """cff-version: 1.2.0
message: If you use this software, please cite it as below.
title: Test Bioinformatics Tool
authors:
  - family-names: Smith
    given-names: Jane
    orcid: https://orcid.org/0000-0000-0000-0001
  - family-names: Doe
    given-names: John
license: MIT
version: 1.0.0
date-released: 2023-01-01
repository-code: https://github.com/example/test-tool
"""

    # Create valid LICENSE file
    license_content = """MIT License

Copyright (c) 2023 Jane Smith, John Doe

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

    # Create valid omnibenchmark.yaml
    omnibenchmark_content = """name: test-bioinformatics-tool
version: 1.0.0
description: A comprehensive test tool for bioinformatics analysis
type: method
stage: preprocessing
entrypoints:
  default: run.sh
input_formats:
  - fastq
  - fasta
output_formats:
  - bam
  - vcf
"""

    # Write files to temporary directory
    (tmp_path / "CITATION.cff").write_text(citation_content)
    (tmp_path / "LICENSE").write_text(license_content)
    (tmp_path / "omnibenchmark.yaml").write_text(omnibenchmark_content)
    (tmp_path / "run.sh").write_text("#!/bin/bash\necho 'test'")

    # Test using the comprehensive validation function
    files_present = {
        "CITATION.cff": (tmp_path / "CITATION.cff").exists(),
        "LICENSE": (tmp_path / "LICENSE").exists(),
        "omnibenchmark.yaml": (tmp_path / "omnibenchmark.yaml").exists(),
        "run.sh": (tmp_path / "run.sh").exists(),
    }

    result = validate_module_files(
        module_id="test_module",
        citation_content=citation_content,
        license_content=license_content,
        omnibenchmark_content=omnibenchmark_content,
        files_present=files_present,
        warn_mode=True,  # Use warn mode to collect all validation results
    )

    # All validations should pass
    assert result.is_valid()
    # Print warnings for debugging
    if result.has_warnings():
        print(f"Warnings found: {result.warnings}")
    assert not result.has_warnings()

    # Test individual components too
    ctx = create_validation_context("test_module", warn_mode=True)
    citation_result, citation_data = validate_citation_cff_content(
        citation_content, ctx
    )

    assert citation_result.is_valid()
    assert citation_data is not None
    assert citation_data["title"] == "Test Bioinformatics Tool"
    assert len(citation_data["authors"]) == 2
    assert citation_data["license"] == "MIT"


@pytest.mark.short
def test_validate_module_with_issues(tmp_path):
    """Test validation of a module with various issues using temporary files."""
    # Create CITATION.cff with multiple issues
    citation_content = """cff-version: 2.0.0
# Missing required 'message' field
title: Test Tool
authors:
  - given-names: John  # Missing family-names
  - family-names: Smith
    # Missing given-names (warning)
license: INVALID_LICENSE
"""

    # Create empty LICENSE file
    license_content = ""

    # Create invalid omnibenchmark.yaml
    omnibenchmark_content = """name: test-tool
invalid: yaml: structure: here
"""

    # Write files to temporary directory
    (tmp_path / "CITATION.cff").write_text(citation_content)
    (tmp_path / "LICENSE").write_text(license_content)
    (tmp_path / "omnibenchmark.yaml").write_text(omnibenchmark_content)

    # Test validation in warn mode to collect all issues
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

    # Should be valid (errors converted to warnings) but have many warnings
    assert result.is_valid()
    assert result.has_warnings()
    assert len(result.warnings) > 0

    # Check for specific issue types
    issue_types = {issue.issue_type for issue in result.issues}
    expected_types = {
        "citation_missing_required_field",
        "citation_unsupported_version",
        "citation_invalid_license",
        "citation_author_missing_family_name",
        "license_file_empty",
        "omnibenchmark_yaml_invalid_yaml",
    }

    # Should have some of the expected issue types
    assert len(issue_types.intersection(expected_types)) > 0

    # Test strict mode - should fail fast on individual validation
    ctx_strict = create_validation_context("test_module", warn_mode=False)
    with pytest.raises(ValidationException):
        validate_citation_cff_content(citation_content, ctx_strict)


@pytest.mark.short
def test_validate_minimal_module(tmp_path):
    """Test validation of a minimal module with only required files."""
    # Create minimal CITATION.cff
    citation_content = """cff-version: 1.2.0
message: Please cite this software
title: Minimal Tool
authors:
  - family-names: Developer
license: Apache-2.0
"""

    # Only create CITATION.cff, no LICENSE or omnibenchmark.yaml
    (tmp_path / "CITATION.cff").write_text(citation_content)

    # Test file structure with missing optional files
    files_present = {
        "CITATION.cff": (tmp_path / "CITATION.cff").exists(),
        "LICENSE": (tmp_path / "LICENSE").exists(),
        "LICENSE.txt": (tmp_path / "LICENSE.txt").exists(),
        "LICENSE.md": (tmp_path / "LICENSE.md").exists(),
        "COPYING": (tmp_path / "COPYING").exists(),
        "COPYING.txt": (tmp_path / "COPYING.txt").exists(),
        "omnibenchmark.yaml": (tmp_path / "omnibenchmark.yaml").exists(),
    }

    result = validate_module_files(
        module_id="test_module",
        citation_content=citation_content,
        files_present=files_present,
        warn_mode=True,
    )

    # omnibenchmark.yaml missing is now an error (recent change), so validation should fail
    assert not result.is_valid()
    assert result.has_warnings()

    # Should have error for missing omnibenchmark.yaml
    error_types = {e.issue_type for e in result.errors}
    assert "omnibenchmark_yaml_missing" in error_types

    # Check for specific warning types
    issue_types = {issue.issue_type for issue in result.warnings}
    # License file warning is suppressed because citation contains license field
    assert "no_license_file" not in issue_types
    # Should have citation author warning since given-names is missing
    assert "citation_author_missing_given_name" in issue_types

    # Test individual citation validation
    ctx = create_validation_context("test_module", warn_mode=True)
    citation_result, citation_data = validate_citation_cff_content(
        citation_content, ctx
    )
    assert citation_result.is_valid()
    assert citation_data["title"] == "Minimal Tool"


@pytest.mark.short
def test_validate_license_consistency_scenarios(tmp_path):
    """Test different license consistency scenarios using temporary files."""
    scenarios = [
        {
            "name": "matching_licenses",
            "citation_license": "MIT",
            "license_content": "MIT License\n\nCopyright (c) 2023...",
            "should_have_warnings": False,
        },
        {
            "name": "mismatched_licenses",
            "citation_license": "MIT",
            "license_content": "Apache License Version 2.0...",
            "should_have_warnings": True,
        },
        {
            "name": "no_license_in_citation",
            "citation_license": None,
            "license_content": "MIT License\n\nCopyright...",
            "should_have_warnings": True,
        },
    ]

    for scenario in scenarios:
        # Test license consistency for each scenario
        ctx = create_validation_context("test_module", warn_mode=True)
        result = validate_license_consistency(
            scenario["citation_license"], scenario["license_content"], ctx
        )

        assert result.is_valid()  # Consistency issues are warnings, not errors

        if scenario["should_have_warnings"]:
            assert (
                result.has_warnings()
            ), f"Scenario '{scenario['name']}' should have warnings"
        else:
            assert (
                not result.has_warnings()
            ), f"Scenario '{scenario['name']}' should not have warnings"


@pytest.mark.short
def test_validate_alternative_license_files(tmp_path):
    """Test validation with alternative LICENSE file names."""
    citation_content = """cff-version: 1.2.0
message: Please cite this software
title: Test Tool
authors:
  - family-names: Developer
license: BSD-3-Clause
"""

    license_content = """BSD 3-Clause License

Copyright (c) 2023, Developer
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
"""

    # Create CITATION.cff and LICENSE.txt (not LICENSE)
    (tmp_path / "CITATION.cff").write_text(citation_content)
    (tmp_path / "LICENSE.txt").write_text(license_content)

    # Test file structure validation
    files_present = {
        "CITATION.cff": True,
        "LICENSE": False,
        "LICENSE.txt": True,  # Alternative license file
        "LICENSE.md": False,
        "COPYING": False,
        "COPYING.txt": False,
        "omnibenchmark.yaml": False,
    }

    ctx = create_validation_context("test_module", warn_mode=True)
    structure_result = validate_file_structure(files_present, ctx)

    # validate_file_structure only checks files, it doesn't error on missing omnibenchmark.yaml
    # (that's done by validate_module_files when validating content)
    assert structure_result.is_valid()

    # Should not warn about LICENSE since LICENSE.txt exists
    issue_types = {issue.issue_type for issue in structure_result.warnings}
    assert "no_license_file" not in issue_types


@pytest.mark.short
def test_validation_with_real_file_reading(tmp_path):
    """Test validation by actually reading files from disk."""
    # Create files
    citation_content = """cff-version: 1.2.0
message: If you use this software, please cite it using these metadata.
title: File Reading Test Tool
authors:
  - family-names: Tester
    given-names: File
license: GPL-3.0-or-later
"""

    (tmp_path / "CITATION.cff").write_text(citation_content)
    (tmp_path / "README.md").write_text("# Test Tool\n\nThis is a test.")

    # Read files and validate
    citation_file = tmp_path / "CITATION.cff"
    readme_file = tmp_path / "README.md"

    assert citation_file.exists()
    assert readme_file.exists()

    # Read and validate CITATION.cff
    actual_content = citation_file.read_text()
    ctx = create_validation_context("test_module", warn_mode=True)
    result, data = validate_citation_cff_content(actual_content, ctx)

    assert result.is_valid()
    assert data["title"] == "File Reading Test Tool"
    assert data["license"] == "GPL-3.0-or-later"

    # Verify file content matches what we wrote
    assert "File Reading Test Tool" in actual_content
    assert "GPL-3.0-or-later" in actual_content


@pytest.mark.short
def test_complex_validation_workflow(tmp_path):
    """Test a complex validation workflow that mimics real usage."""
    # Simulate a validation workflow for a complete module

    # Step 1: Check file structure first
    files_present = {
        "CITATION.cff": (tmp_path / "CITATION.cff").exists(),
        "LICENSE": (tmp_path / "LICENSE").exists(),
        "omnibenchmark.yaml": (tmp_path / "omnibenchmark.yaml").exists(),
    }

    ctx = create_validation_context("test_module", warn_mode=False)

    # Initially should fail - no files exist
    with pytest.raises(ValidationException) as exc_info:
        validate_file_structure(files_present, ctx)

    assert exc_info.value.issues[0].issue_type == "citation_missing"

    # Step 2: Create CITATION.cff
    citation_content = """cff-version: 1.2.0
message: If you use this software, please cite it as below.
title: Complex Validation Tool
authors:
  - family-names: Engineer
    given-names: Expert
license: MIT
"""
    (tmp_path / "CITATION.cff").write_text(citation_content)

    # Re-check structure in warn mode
    files_present["CITATION.cff"] = True
    ctx_warn = create_validation_context("test_module", warn_mode=True)
    structure_result = validate_file_structure(files_present, ctx_warn)
    assert structure_result.is_valid()  # Now valid, but with warnings
    assert structure_result.has_warnings()

    # Step 3: Validate CITATION.cff content
    ctx_cite = create_validation_context("test_module", warn_mode=True)
    citation_result, citation_data = validate_citation_cff_content(
        citation_content, ctx_cite
    )
    assert citation_result.is_valid()

    # Step 4: Add LICENSE file
    license_content = "MIT License\n\nCopyright (c) 2023 Expert Engineer\n..."
    (tmp_path / "LICENSE").write_text(license_content)

    # Step 5: Test license consistency
    ctx_license = create_validation_context("test_module", warn_mode=True)
    consistency_result = validate_license_consistency(
        citation_data["license"], license_content, ctx_license
    )
    assert consistency_result.is_valid()

    # Step 6: Final validation - everything should be good
    files_present["LICENSE"] = True
    final_result = validate_module_files(
        module_id="test_module",
        citation_content=citation_content,
        license_content=license_content,
        files_present=files_present,
        warn_mode=True,
    )

    # omnibenchmark.yaml missing is now an error, so validation should fail
    assert not final_result.is_valid()

    # Should have error for missing omnibenchmark.yaml
    error_types = {e.issue_type for e in final_result.errors}
    assert "omnibenchmark_yaml_missing" in error_types

    # Should not warn about license
    issue_types = {issue.issue_type for issue in final_result.warnings}
    assert "no_license_file" not in issue_types


@pytest.mark.short
def test_strategy_pattern_integration():
    """Test that the strategy pattern works correctly for different modes."""
    invalid_citation = """cff-version: 1.2.0
# Missing: message, title, authors
license: INVALID-LICENSE
"""

    # Test strict mode - should fail fast
    ctx_strict = create_validation_context("test_module", warn_mode=False)
    with pytest.raises(ValidationException) as exc_info:
        validate_citation_cff_content(invalid_citation, ctx_strict)

    assert len(exc_info.value.issues) == 1
    assert exc_info.value.issues[0].is_error

    # Test warn mode - should collect all issues as warnings
    ctx_warn = create_validation_context("test_module", warn_mode=True)
    result, data = validate_citation_cff_content(invalid_citation, ctx_warn)

    assert result.is_valid()  # Valid because errors converted to warnings
    assert len(ctx_warn.result.warnings) > 0
    assert all(issue.is_warning for issue in ctx_warn.result.issues)


@pytest.mark.short
def test_deduplication_and_aggregation():
    """Test that the validation system properly deduplicates and aggregates issues."""
    # Create validation context
    ctx = create_validation_context("test_module", warn_mode=True)

    # Add the same issue multiple times
    for i in range(3):
        ctx.add_issue(
            issue_type="test_duplicate",
            path="test.txt",
            message="This is a duplicate issue",
        )

    # Should only have one issue due to deduplication
    assert len(ctx.result.issues) == 1

    # Add issues for different modules

    result = ValidationResult()

    # Module 1 issues
    result.add_issue(
        ValidationIssue(
            issue_type="citation_missing",
            severity=ValidationSeverity.ERROR,
            path="module1/CITATION.cff",
            module_id="module1",
            message="CITATION.cff not found",
        )
    )

    # Module 2 issues
    result.add_issue(
        ValidationIssue(
            issue_type="license_missing",
            severity=ValidationSeverity.WARNING,
            path="module2/LICENSE",
            module_id="module2",
            message="LICENSE file not found",
        )
    )

    # Test aggregation by module
    by_module = result.get_issues_by_module()
    assert len(by_module) == 2
    assert "module1" in by_module
    assert "module2" in by_module
    assert len(by_module["module1"]) == 1
    assert len(by_module["module2"]) == 1
