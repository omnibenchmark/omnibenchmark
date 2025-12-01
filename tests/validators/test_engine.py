"""Unit tests for the artifact validation engine."""

import json
from dataclasses import dataclass
from pathlib import Path

import pytest
import yaml
from returns.result import Success, Failure

from omnibenchmark.validators import (
    run_engine,
    run_validation_for_module,
    apply_rule,
    load_validation_rules,
    ValidationResult,
    ModuleValidationResult,
    ArtifactNotFoundError,
    VALIDATOR_REGISTRY,
)


@dataclass(frozen=True, slots=True)
class MockModule:
    """Test implementation of ValidatableModule protocol."""

    name: str
    path: Path
    artifact_paths: dict[str, Path]


@dataclass(frozen=True, slots=True)
class MockBenchmark:
    """Test implementation of ValidatableBenchmark protocol."""

    path: Path
    modules: list[MockModule]

    def __iter__(self):
        yield from self.modules


pytestmark = pytest.mark.short


@pytest.fixture
def temp_module_dir(tmp_path):
    """Create a temporary module directory with test artifacts."""
    module_path = tmp_path / "test_module"
    module_path.mkdir()

    # Create test artifacts
    (module_path / "output.txt").write_text("Hello World\nTest Content")
    (module_path / "results.csv").write_text("id,score,value\n1,10.5,foo\n2,20.3,bar\n")
    (module_path / "metadata.json").write_text(
        json.dumps({"timestamp": "2025-01-01", "version": "1.0"})
    )
    (module_path / "data_dir").mkdir()

    return module_path


@pytest.fixture
def validation_yaml_path(temp_module_dir):
    """Create a validation.yaml file."""
    validation_content = {
        "output.txt": [
            {"is_file": True},
            {"matches": "Hello"},
        ],
        "results.csv": [
            {"is_file": True},
            {"has_column": "score"},
            {"not_na": "score"},
        ],
        "metadata.json": [
            {"is_file": True},
            {"has_key": "timestamp"},
        ],
        "data_dir": [
            {"is_dir": True},
        ],
    }

    validation_file = temp_module_dir / "validation.yaml"
    with open(validation_file, "w") as f:
        yaml.dump(validation_content, f)

    return validation_file


class TestValidatorRegistry:
    """Test the validator registry."""

    def test_registry_has_all_validators(self):
        """Ensure all expected validators are registered."""
        expected_validators = {
            "is_file",
            "is_dir",
            "matches",
            "checksum",
            "has_column",
            "not_na",
            "has_key",
        }
        assert set(VALIDATOR_REGISTRY.keys()) == expected_validators


class TestIsFileValidator:
    """Test the is_file validator."""

    def test_is_file_success(self, temp_module_dir):
        """Test is_file passes for actual files."""
        rule = {"is_file": True}
        result = apply_rule(
            "output.txt", temp_module_dir / "output.txt", rule, temp_module_dir
        )
        assert result.success is True
        assert "file" in result.message.lower()

    def test_is_file_fails_for_directory(self, temp_module_dir):
        """Test is_file fails for directories."""
        rule = {"is_file": True}
        result = apply_rule(
            "data_dir", temp_module_dir / "data_dir", rule, temp_module_dir
        )
        assert result.success is False
        assert "not a file" in result.message.lower()

    def test_is_file_fails_for_missing(self, temp_module_dir):
        """Test is_file fails for missing files."""
        rule = {"is_file": True}
        result = apply_rule(
            "missing.txt", temp_module_dir / "missing.txt", rule, temp_module_dir
        )
        assert result.success is False
        assert "does not exist" in result.message.lower()


class TestIsDirValidator:
    """Test the is_dir validator."""

    def test_is_dir_success(self, temp_module_dir):
        """Test is_dir passes for actual directories."""
        rule = {"is_dir": True}
        result = apply_rule(
            "data_dir", temp_module_dir / "data_dir", rule, temp_module_dir
        )
        assert result.success is True
        assert "directory" in result.message.lower()

    def test_is_dir_fails_for_file(self, temp_module_dir):
        """Test is_dir fails for files."""
        rule = {"is_dir": True}
        result = apply_rule(
            "output.txt", temp_module_dir / "output.txt", rule, temp_module_dir
        )
        assert result.success is False
        assert "not a directory" in result.message.lower()


class TestMatchesValidator:
    """Test the regex matches validator."""

    def test_matches_success(self, temp_module_dir):
        """Test matches passes when pattern found."""
        rule = {"matches": "Hello"}
        result = apply_rule(
            "output.txt", temp_module_dir / "output.txt", rule, temp_module_dir
        )
        assert result.success is True

    def test_matches_with_regex(self, temp_module_dir):
        """Test matches works with regex patterns."""
        rule = {"matches": r"Test\s+Content"}
        result = apply_rule(
            "output.txt", temp_module_dir / "output.txt", rule, temp_module_dir
        )
        assert result.success is True

    def test_matches_failure(self, temp_module_dir):
        """Test matches fails when pattern not found."""
        rule = {"matches": "NotPresent"}
        result = apply_rule(
            "output.txt", temp_module_dir / "output.txt", rule, temp_module_dir
        )
        assert result.success is False
        assert "does not match" in result.message.lower()


class TestChecksumValidator:
    """Test the checksum validator."""

    def test_checksum_md5_success(self, temp_module_dir):
        """Test checksum passes with correct MD5 hash."""
        import hashlib

        content = (temp_module_dir / "output.txt").read_bytes()
        expected_hash = hashlib.md5(content).hexdigest()

        rule = {"checksum": f"md5:{expected_hash}"}
        result = apply_rule(
            "output.txt", temp_module_dir / "output.txt", rule, temp_module_dir
        )
        assert result.success is True
        assert "checksum" in result.message.lower()

    def test_checksum_sha256_success(self, temp_module_dir):
        """Test checksum passes with correct SHA256 hash."""
        import hashlib

        content = (temp_module_dir / "output.txt").read_bytes()
        expected_hash = hashlib.sha256(content).hexdigest()

        rule = {"checksum": f"sha256:{expected_hash}"}
        result = apply_rule(
            "output.txt", temp_module_dir / "output.txt", rule, temp_module_dir
        )
        assert result.success is True

    def test_checksum_failure(self, temp_module_dir):
        """Test checksum fails with wrong hash."""
        rule = {"checksum": "md5:wronghash123"}
        result = apply_rule(
            "output.txt", temp_module_dir / "output.txt", rule, temp_module_dir
        )
        assert result.success is False
        assert "mismatch" in result.message.lower()


class TestHasColumnValidator:
    """Test the CSV has_column validator."""

    def test_has_column_success(self, temp_module_dir):
        """Test has_column passes when column exists."""
        rule = {"has_column": "score"}
        result = apply_rule(
            "results.csv", temp_module_dir / "results.csv", rule, temp_module_dir
        )
        assert result.success is True
        assert "found" in result.message.lower()

    def test_has_column_failure(self, temp_module_dir):
        """Test has_column fails when column missing."""
        rule = {"has_column": "missing_column"}
        result = apply_rule(
            "results.csv", temp_module_dir / "results.csv", rule, temp_module_dir
        )
        assert result.success is False
        assert "not found" in result.message.lower()


class TestNotNAValidator:
    """Test the CSV not_na validator."""

    def test_not_na_success(self, temp_module_dir):
        """Test not_na passes when column has no NA values."""
        rule = {"not_na": "score"}
        result = apply_rule(
            "results.csv", temp_module_dir / "results.csv", rule, temp_module_dir
        )
        assert result.success is True
        assert "no n/a" in result.message.lower()

    def test_not_na_failure_empty(self, temp_module_dir):
        """Test not_na fails when column has empty values."""
        # Create CSV with empty value
        (temp_module_dir / "bad.csv").write_text("id,score\n1,10\n2,\n")

        rule = {"not_na": "score"}
        result = apply_rule(
            "bad.csv", temp_module_dir / "bad.csv", rule, temp_module_dir
        )
        assert result.success is False
        assert "contains n/a" in result.message.lower()

    def test_not_na_failure_na_string(self, temp_module_dir):
        """Test not_na fails when column has 'NA' string."""
        # Create CSV with NA string
        (temp_module_dir / "bad.csv").write_text("id,score\n1,10\n2,NA\n")

        rule = {"not_na": "score"}
        result = apply_rule(
            "bad.csv", temp_module_dir / "bad.csv", rule, temp_module_dir
        )
        assert result.success is False

    def test_not_na_missing_column(self, temp_module_dir):
        """Test not_na fails gracefully when column doesn't exist."""
        rule = {"not_na": "missing_column"}
        result = apply_rule(
            "results.csv", temp_module_dir / "results.csv", rule, temp_module_dir
        )
        assert result.success is False
        assert "not found" in result.message.lower()


class TestHasKeyValidator:
    """Test the JSON has_key validator."""

    def test_has_key_success(self, temp_module_dir):
        """Test has_key passes when key exists."""
        rule = {"has_key": "timestamp"}
        result = apply_rule(
            "metadata.json", temp_module_dir / "metadata.json", rule, temp_module_dir
        )
        assert result.success is True
        assert "found" in result.message.lower()

    def test_has_key_failure(self, temp_module_dir):
        """Test has_key fails when key missing."""
        rule = {"has_key": "missing_key"}
        result = apply_rule(
            "metadata.json", temp_module_dir / "metadata.json", rule, temp_module_dir
        )
        assert result.success is False
        assert "not found" in result.message.lower()

    def test_has_key_invalid_json(self, temp_module_dir):
        """Test has_key handles non-object JSON."""
        # Create JSON array instead of object
        (temp_module_dir / "array.json").write_text("[1, 2, 3]")

        rule = {"has_key": "foo"}
        result = apply_rule(
            "array.json", temp_module_dir / "array.json", rule, temp_module_dir
        )
        assert result.success is False
        assert "not an object" in result.message.lower()


class TestLoadValidationRules:
    """Test loading validation rules from YAML."""

    def test_load_rules_success(self, temp_module_dir, validation_yaml_path):
        """Test loading validation rules successfully."""
        result = load_validation_rules(temp_module_dir)

        assert isinstance(result, Success)
        rules = result.unwrap()
        assert "output.txt" in rules
        assert "results.csv" in rules
        assert len(rules["output.txt"]) == 2

    def test_load_rules_missing_file(self, temp_module_dir):
        """Test loading rules fails when validation.yaml missing."""
        result = load_validation_rules(temp_module_dir)

        assert isinstance(result, Failure)
        error = result.failure()
        assert isinstance(error, ArtifactNotFoundError)


class TestModuleValidation:
    """Test validation at the module level."""

    def test_run_validation_for_module_success(
        self, temp_module_dir, validation_yaml_path
    ):
        """Test running validation for a module with all passing rules."""
        artifact_paths = {
            "output.txt": temp_module_dir / "output.txt",
            "results.csv": temp_module_dir / "results.csv",
            "metadata.json": temp_module_dir / "metadata.json",
            "data_dir": temp_module_dir / "data_dir",
        }

        module = MockModule(
            name="test_module", path=temp_module_dir, artifact_paths=artifact_paths
        )

        result = run_validation_for_module(module)

        assert result.module_name == "test_module"
        assert result.error is None
        assert result.all_successful is True
        assert len(result.results) > 0

    def test_run_validation_missing_artifact(
        self, temp_module_dir, validation_yaml_path
    ):
        """Test validation fails when artifact path not provided."""
        artifact_paths = {
            "output.txt": temp_module_dir / "output.txt",
            # Missing other artifacts
        }

        module = MockModule(
            name="test_module", path=temp_module_dir, artifact_paths=artifact_paths
        )

        result = run_validation_for_module(module)

        assert result.all_successful is False
        # Should have failures for missing artifacts
        failures = [r for r in result.results if not r.success]
        assert len(failures) > 0

    def test_run_validation_no_yaml(self, temp_module_dir):
        """Test validation when no validation.yaml exists."""
        module = MockModule(name="test_module", path=temp_module_dir, artifact_paths={})

        result = run_validation_for_module(module)

        assert result.error is not None
        assert "not found" in result.error.lower()


class TestEngineOrchestration:
    """Test the full validation engine."""

    def test_run_engine(self, temp_module_dir, validation_yaml_path):
        """Test running the engine on a benchmark with multiple modules."""
        # Create another module
        module2_path = temp_module_dir.parent / "module2"
        module2_path.mkdir()
        (module2_path / "output.txt").write_text("Module 2 output")

        validation_content = {"output.txt": [{"is_file": True}]}
        with open(module2_path / "validation.yaml", "w") as f:
            yaml.dump(validation_content, f)

        modules = [
            MockModule(
                name="module1",
                path=temp_module_dir,
                artifact_paths={
                    "output.txt": temp_module_dir / "output.txt",
                    "results.csv": temp_module_dir / "results.csv",
                    "metadata.json": temp_module_dir / "metadata.json",
                    "data_dir": temp_module_dir / "data_dir",
                },
            ),
            MockModule(
                name="module2",
                path=module2_path,
                artifact_paths={
                    "output.txt": module2_path / "output.txt",
                },
            ),
        ]

        benchmark = MockBenchmark(path=temp_module_dir.parent, modules=modules)

        results = run_engine(benchmark)

        assert len(results) == 2
        assert results[0].module_name == "module1"
        assert results[1].module_name == "module2"
        assert all(r.all_successful for r in results)


class TestValidationResult:
    """Test ValidationResult data class."""

    def test_validation_result_creation(self):
        """Test creating a ValidationResult."""
        rule = {"is_file": True}
        result = ValidationResult(
            artifact_name="test.txt", rule=rule, success=True, message="Test passed"
        )

        assert result.artifact_name == "test.txt"
        assert result.success is True
        assert result.message == "Test passed"


class TestModuleValidationResult:
    """Test ModuleValidationResult data class."""

    def test_all_successful_true(self):
        """Test all_successful returns True when all pass."""
        results = [
            ValidationResult("a", {"is_file": True}, True, "ok"),
            ValidationResult("b", {"is_file": True}, True, "ok"),
        ]

        module_result = ModuleValidationResult(module_name="test", results=results)

        assert module_result.all_successful is True

    def test_all_successful_false(self):
        """Test all_successful returns False when any fail."""
        results = [
            ValidationResult("a", {"is_file": True}, True, "ok"),
            ValidationResult("b", {"is_file": True}, False, "fail"),
        ]

        module_result = ModuleValidationResult(module_name="test", results=results)

        assert module_result.all_successful is False

    def test_all_successful_with_error(self):
        """Test all_successful returns False when error present."""
        module_result = ModuleValidationResult(
            module_name="test", results=[], error="Something went wrong"
        )

        assert module_result.all_successful is False
