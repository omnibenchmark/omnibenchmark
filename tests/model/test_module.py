import pytest

from pydantic import ValidationError

from omnibenchmark.model.module import ModuleMetadata, DerivedSoftware

VALID_YAML = """
name: "Test Module"
description: "A test module for omnibenchmark"
author: "Test Author"
derives_from:
  - description: "Source package"
    author: "Original Author"
    url: "https://example.com/source"
    doi: "10.1234/example.doi"
license: "MIT"
"""


@pytest.mark.short
def test_valid_module_from_yaml():
    """Test creating a ModuleMetadata from valid YAML"""
    module = ModuleMetadata.from_yaml(VALID_YAML)
    assert module.name == "Test Module"
    assert module.description == "A test module for omnibenchmark"
    assert module.author == "Test Author"
    assert module.license == "MIT"
    assert len(module.derives_from) == 1
    assert module.derives_from[0].doi == "10.1234/example.doi"


@pytest.mark.short
def test_license_validation():
    """Test that license validation works with known licenses"""
    valid_module = ModuleMetadata(
        name="Test",
        description="Test",
        author="Test",
        derives_from=[],
        license="MIT",  # Known license
    )
    valid_module.validate_license()  # Should not raise

    invalid_module = ModuleMetadata(
        name="Test",
        description="Test",
        author="Test",
        derives_from=[],
        license="NOT_A_LICENSE",  # Unknown license
    )
    with pytest.raises(ValueError, match="not in recognized open source licenses"):
        invalid_module.validate_license()


@pytest.mark.short
def test_derived_software_validation():
    """Test validation of DerivedSoftware fields"""
    # Valid case
    valid_derived = DerivedSoftware(
        description="Test",
        author="Test",
        url="https://example.com",
        doi="10.1234/example.doi",
    )
    assert valid_derived.doi == "10.1234/example.doi"

    # Invalid URL
    with pytest.raises(ValidationError):
        DerivedSoftware(
            description="Test",
            author="Test",
            url="not-a-url",
            doi="10.1234/example.doi",
        )


@pytest.mark.short
def test_missing_required_fields():
    """Test that required fields are enforced"""
    incomplete_yaml = """
    name: "Incomplete Module"
    description: "Missing required fields"
    """
    with pytest.raises(ValidationError):
        ModuleMetadata.from_yaml(incomplete_yaml)


@pytest.mark.short
def test_optional_fields():
    """Test that optional fields work as expected"""
    yaml_with_optionals = """
    name: "Module with optionals"
    description: "Test optional fields"
    author: "Test Author"
    derives_from:
      - description: "Source with easyconfig"
        author: "Original Author"
        url: "https://example.com/source"
        easyconfig: "file:///path/to/config"
    license: "Apache-2.0"
    """
    module = ModuleMetadata.from_yaml(yaml_with_optionals)
    assert module.derives_from[0].easyconfig == "file:///path/to/config"
    assert module.derives_from[0].doi is None
