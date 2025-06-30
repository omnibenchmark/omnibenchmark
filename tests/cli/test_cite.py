import pytest
import tempfile
import shutil
from pathlib import Path
from unittest.mock import Mock, patch

from omnibenchmark.benchmark.cite import CitationExtractionError, create_cite_issue
from tests.cli.cli_setup import OmniCLISetup


@pytest.fixture
def cli_setup():
    """Set up CLI testing environment."""
    return OmniCLISetup()


@pytest.fixture
def temp_benchmark_dir():
    """Create a temporary benchmark directory for testing."""
    temp_dir = Path(tempfile.mkdtemp())
    yield temp_dir
    shutil.rmtree(temp_dir)


@pytest.fixture
def mock_benchmark_config(temp_benchmark_dir):
    """Create a mock benchmark configuration file."""
    config_content = """
benchmark:
  name: test_benchmark
  modules:
    - name: test_module
      repository: https://github.com/test/repo.git
"""
    config_file = temp_benchmark_dir / "config.yaml"
    config_file.write_text(config_content)
    return config_file


class TestCiteCLI:
    """Test suite for the CLI cite command."""

    @pytest.mark.short
    def test_cite_help(self, cli_setup):
        """Test that cite command shows help information."""
        result = cli_setup.call(["info", "cite", "--help"])

        assert result.returncode == 0
        assert "Extract citation metadata" in result.stdout
        assert "--format" in result.stdout
        assert "--out" in result.stdout
        assert "--out" in result.stdout
        assert "--warn" in result.stdout

    @pytest.mark.short
    def test_cite_missing_benchmark(self, cli_setup, tmp_path):
        """Test cite command with missing benchmark file."""
        # Change to a directory without benchmark config
        result = cli_setup.call(["info", "cite"], cwd=str(tmp_path))

        assert result.returncode != 0
        # Should fail due to missing benchmark configuration

    @pytest.mark.short
    @patch("omnibenchmark.cli.benchmark.extract_citation_metadata")
    @patch("omnibenchmark.cli.utils.validation.validate_benchmark")
    def test_cite_success_json_format(
        self,
        mock_validate,
        mock_extract,
        cli_setup,
        temp_benchmark_dir,
        mock_benchmark_config,
    ):
        """Test successful citation extraction with JSON format."""
        # Mock the benchmark validation
        mock_benchmark = Mock()
        mock_validate.return_value = mock_benchmark

        # Mock successful citation extraction
        citation_data = {
            "test_module": {
                "module_id": "test_module",
                "repository_url": "https://github.com/test/repo.git",
                "commit_hash": "abc123",
                "citation_data": {
                    "title": "Test Software",
                    "authors": [{"family-names": "Doe", "given-names": "John"}],
                },
                "citation_file_found": True,
                "local_repo_exists": True,
            }
        }
        mock_extract.return_value = citation_data

        result = cli_setup.call(
            ["info", "cite", "--format", "json"], cwd=str(temp_benchmark_dir)
        )

        assert result.returncode == 0
        assert "test_module" in result.stdout
        assert "Test Software" in result.stdout

        # Verify the extraction was called with correct parameters
        mock_extract.assert_called_once_with(
            mock_benchmark, strict=True, warn_mode=False
        )

    @pytest.mark.short
    @patch("omnibenchmark.cli.benchmark.extract_citation_metadata")
    @patch("omnibenchmark.cli.utils.validation.validate_benchmark")
    def test_cite_success_bibtex_format(
        self,
        mock_validate,
        mock_extract,
        cli_setup,
        temp_benchmark_dir,
        mock_benchmark_config,
    ):
        """Test successful citation extraction with BibTeX format."""
        # Mock the benchmark validation
        mock_benchmark = Mock()
        mock_validate.return_value = mock_benchmark

        # Mock successful citation extraction
        citation_data = {
            "test_module": {
                "module_id": "test_module",
                "repository_url": "https://github.com/test/repo.git",
                "commit_hash": "abc123",
                "citation_data": {
                    "title": "Test Software",
                    "authors": [{"family-names": "Doe", "given-names": "John"}],
                },
                "citation_file_found": True,
                "local_repo_exists": True,
            }
        }
        mock_extract.return_value = citation_data

        result = cli_setup.call(
            ["info", "cite", "--format", "bibtex"], cwd=str(temp_benchmark_dir)
        )

        assert result.returncode == 0
        # BibTeX format should contain typical BibTeX structure
        assert "@" in result.stdout  # BibTeX entries start with @

    @pytest.mark.short
    @patch("omnibenchmark.cli.benchmark.extract_citation_metadata")
    @patch("omnibenchmark.cli.utils.validation.validate_benchmark")
    def test_cite_output_to_file(
        self,
        mock_validate,
        mock_extract,
        cli_setup,
        temp_benchmark_dir,
        mock_benchmark_config,
    ):
        """Test citation extraction with output to file."""
        # Mock the benchmark validation
        mock_benchmark = Mock()
        mock_validate.return_value = mock_benchmark

        # Mock successful citation extraction
        citation_data = {
            "test_module": {
                "module_id": "test_module",
                "repository_url": "https://github.com/test/repo.git",
                "commit_hash": "abc123",
                "citation_data": {
                    "title": "Test Software",
                    "authors": [{"family-names": "Doe", "given-names": "John"}],
                },
                "citation_file_found": True,
                "local_repo_exists": True,
            }
        }
        mock_extract.return_value = citation_data

        output_file = temp_benchmark_dir / "citations.json"

        result = cli_setup.call(
            ["info", "cite", "--format", "json", "--out", str(output_file)],
            cwd=str(temp_benchmark_dir),
        )

        assert result.returncode == 0
        assert output_file.exists()

        # Verify file content
        content = output_file.read_text()
        assert "test_module" in content
        assert "Test Software" in content

    @pytest.mark.short
    @patch("omnibenchmark.cli.benchmark.extract_citation_metadata")
    @patch("omnibenchmark.cli.utils.validation.validate_benchmark")
    def test_cite_citation_extraction_error(
        self,
        mock_validate,
        mock_extract,
        cli_setup,
        temp_benchmark_dir,
        mock_benchmark_config,
    ):
        """Test handling of CitationExtractionError in strict mode."""
        # Mock the benchmark validation
        mock_benchmark = Mock()
        mock_validate.return_value = mock_benchmark

        # Mock citation extraction error
        issues = [
            create_cite_issue(
                "missing_repository_info",
                "test_module",
                "Missing repository information",
            ),
            create_cite_issue(
                "citation_missing", "test_module", "Citation file not found"
            ),
        ]

        mock_extract.side_effect = CitationExtractionError(
            "Citation extraction failed for module: test_module",
            ["test_module"],
            issues,
        )

        result = cli_setup.call(["info", "cite"], cwd=str(temp_benchmark_dir))

        assert result.returncode == 1
        assert "Citation extraction failed" in result.stderr
        assert "Failed modules: test_module" in result.stderr
        assert "Missing repository information" in result.stderr
        assert "Citation file not found" in result.stderr

    @pytest.mark.short
    @patch("omnibenchmark.cli.benchmark.extract_citation_metadata")
    @patch("omnibenchmark.cli.utils.validation.validate_benchmark")
    def test_cite_runtime_warning(
        self,
        mock_validate,
        mock_extract,
        cli_setup,
        temp_benchmark_dir,
        mock_benchmark_config,
    ):
        """Test handling of RuntimeWarning in warn mode."""
        # Mock the benchmark validation
        mock_benchmark = Mock()
        mock_validate.return_value = mock_benchmark

        # Mock runtime warning
        mock_extract.side_effect = RuntimeWarning("No local repositories found")

        result = cli_setup.call(["info", "cite", "--warn"], cwd=str(temp_benchmark_dir))

        assert result.returncode == 0
        # Warning should be logged but not cause failure

    @pytest.mark.short
    @patch("omnibenchmark.cli.benchmark.extract_citation_metadata")
    @patch("omnibenchmark.cli.utils.validation.validate_benchmark")
    def test_cite_invalid_format(
        self,
        mock_validate,
        mock_extract,
        cli_setup,
        temp_benchmark_dir,
        mock_benchmark_config,
    ):
        """Test handling of invalid output format."""
        # Mock the benchmark validation
        mock_benchmark = Mock()
        mock_validate.return_value = mock_benchmark

        # Mock successful citation extraction
        citation_data = {"test_module": {}}
        mock_extract.return_value = citation_data

        # Mock format_output raising ValueError for invalid format
        with patch("omnibenchmark.cli.benchmark.format_output") as mock_format:
            mock_format.side_effect = ValueError("Unsupported format: invalid")

            result = cli_setup.call(
                ["info", "cite", "--format", "invalid"], cwd=str(temp_benchmark_dir)
            )

            assert result.returncode == 1
            assert "Unsupported format" in result.stderr

    @pytest.mark.short
    @patch("omnibenchmark.cli.benchmark.extract_citation_metadata")
    @patch("omnibenchmark.cli.utils.validation.validate_benchmark")
    def test_cite_file_write_error(
        self,
        mock_validate,
        mock_extract,
        cli_setup,
        temp_benchmark_dir,
        mock_benchmark_config,
    ):
        """Test handling of file write errors."""
        # Mock the benchmark validation
        mock_benchmark = Mock()
        mock_validate.return_value = mock_benchmark

        # Mock successful citation extraction
        citation_data = {"test_module": {}}
        mock_extract.return_value = citation_data

        # Try to write to an invalid path
        invalid_path = "/invalid/path/citations.json"

        result = cli_setup.call(
            ["info", "cite", "--format", "json", "--out", invalid_path],
            cwd=str(temp_benchmark_dir),
        )

        assert result.returncode == 1
        assert "Failed to write output file" in result.stderr

    @pytest.mark.short
    @patch("omnibenchmark.cli.benchmark.extract_citation_metadata")
    @patch("omnibenchmark.cli.utils.validation.validate_benchmark")
    def test_cite_strict_mode_parameter(
        self,
        mock_validate,
        mock_extract,
        cli_setup,
        temp_benchmark_dir,
        mock_benchmark_config,
    ):
        """Test that strict mode parameter is passed correctly."""
        # Mock the benchmark validation
        mock_benchmark = Mock()
        mock_validate.return_value = mock_benchmark

        # Mock successful citation extraction
        citation_data = {"test_module": {}}
        mock_extract.return_value = citation_data

        # Test with --warn flag (which is the opposite of strict)
        result = cli_setup.call(["info", "cite", "--warn"], cwd=str(temp_benchmark_dir))

        assert result.returncode == 0
        # Verify extraction was called with warn_mode=True
        mock_extract.assert_called_with(mock_benchmark, strict=True, warn_mode=True)

    @pytest.mark.short
    @patch("omnibenchmark.cli.benchmark.extract_citation_metadata")
    @patch("omnibenchmark.cli.utils.validation.validate_benchmark")
    def test_cite_warn_mode_parameter(
        self,
        mock_validate,
        mock_extract,
        cli_setup,
        temp_benchmark_dir,
        mock_benchmark_config,
    ):
        """Test that warn mode parameter is passed correctly."""
        # Mock the benchmark validation
        mock_benchmark = Mock()
        mock_validate.return_value = mock_benchmark

        # Mock successful citation extraction
        citation_data = {"test_module": {}}
        mock_extract.return_value = citation_data

        # Test with --warn flag
        result = cli_setup.call(["info", "cite", "--warn"], cwd=str(temp_benchmark_dir))

        assert result.returncode == 0
        # Verify extraction was called with warn_mode=True
        mock_extract.assert_called_with(mock_benchmark, strict=True, warn_mode=True)

    @pytest.mark.short
    @patch("omnibenchmark.cli.benchmark.extract_citation_metadata")
    @patch("omnibenchmark.cli.utils.validation.validate_benchmark")
    def test_cite_multiple_module_errors(
        self,
        mock_validate,
        mock_extract,
        cli_setup,
        temp_benchmark_dir,
        mock_benchmark_config,
    ):
        """Test handling of multiple module extraction errors."""
        # Mock the benchmark validation
        mock_benchmark = Mock()
        mock_validate.return_value = mock_benchmark

        # Mock citation extraction error with multiple modules
        issues = [
            create_cite_issue(
                "missing_repository_info",
                "module1",
                "Missing repository information for module1",
            ),
            create_cite_issue(
                "citation_missing", "module2", "Citation file not found for module2"
            ),
            create_cite_issue(
                "clone_failed",
                "module3",
                "Failed to clone repository for module3",
                repo_url="https://github.com/test/repo3.git",
            ),
        ]

        mock_extract.side_effect = CitationExtractionError(
            "Citation extraction failed for modules: module1, module2, module3",
            ["module1", "module2", "module3"],
            issues,
        )

        result = cli_setup.call(["info", "cite"], cwd=str(temp_benchmark_dir))

        assert result.returncode == 1
        assert "Failed modules: module1, module2, module3" in result.stderr
        assert "Missing repository information for module1" in result.stderr
        assert "Citation file not found for module2" in result.stderr
        assert "Failed to clone repository for module3" in result.stderr

    @pytest.mark.short
    @patch("omnibenchmark.cli.benchmark.extract_citation_metadata")
    @patch("omnibenchmark.cli.utils.validation.validate_benchmark")
    def test_cite_empty_citation_data(
        self,
        mock_validate,
        mock_extract,
        cli_setup,
        temp_benchmark_dir,
        mock_benchmark_config,
    ):
        """Test handling of empty citation data."""
        # Mock the benchmark validation
        mock_benchmark = Mock()
        mock_validate.return_value = mock_benchmark

        # Mock empty citation extraction result
        mock_extract.return_value = {}

        result = cli_setup.call(
            ["info", "cite", "--format", "json"], cwd=str(temp_benchmark_dir)
        )

        assert result.returncode == 0
        assert "{}" in result.stdout or "[]" in result.stdout  # Empty JSON object/array

    @pytest.mark.short
    def test_cite_integration_error_attribute_fix(self, cli_setup, temp_benchmark_dir):
        """Test that the CitationExtractionError attribute fix works in integration."""
        # This test specifically targets the bug that was fixed:
        # AttributeError: 'CitationExtractionError' object has no attribute 'errors'

        # Create a minimal config that will likely trigger citation extraction errors
        config_content = """
benchmark:
  name: test_benchmark
  modules:
    - name: nonexistent_module
      repository: https://github.com/nonexistent/repo.git
"""
        config_file = temp_benchmark_dir / "config.yaml"
        config_file.write_text(config_content)

        # Run the command - this should not crash with AttributeError
        result = cli_setup.call(["info", "cite"], cwd=str(temp_benchmark_dir))

        # The command should fail (likely due to missing repos), but NOT with AttributeError
        # The specific error we're testing against should not appear in stderr
        assert "AttributeError" not in result.stderr
        assert (
            "'CitationExtractionError' object has no attribute 'errors'"
            not in result.stderr
        )

        # Instead, we should see proper error handling
        if result.returncode != 0:
            # Should show proper error messages, not attribute errors
            assert (
                "Failed modules:" in result.stderr
                or "Citation extraction" in result.stderr
            )
