"""
Unit tests for entrypoint reading in execution.py.

Tests the _read_entrypoint function which supports both:
- New-style omnibenchmark.yaml (preferred)
- Old-style config.cfg (deprecated, with FutureWarning)
"""

import pytest
import warnings
from omnibenchmark.workflow.snakemake.scripts.execution import _read_entrypoint


@pytest.mark.short
def test_read_entrypoint_from_omnibenchmark_yaml(tmp_path):
    """Test reading entrypoint from omnibenchmark.yaml."""
    module_dir = tmp_path / "test_module"
    module_dir.mkdir()

    # Create omnibenchmark.yaml
    yaml_file = module_dir / "omnibenchmark.yaml"
    yaml_file.write_text("""# OmniBenchmark Module Configuration

entrypoints:
  default: run.py
""")

    entrypoint = _read_entrypoint(module_dir, "test_module")
    assert entrypoint == "run.py"


@pytest.mark.short
def test_read_entrypoint_from_omnibenchmark_yaml_shell_script(tmp_path):
    """Test reading shell script entrypoint from omnibenchmark.yaml."""
    module_dir = tmp_path / "test_module"
    module_dir.mkdir()

    yaml_file = module_dir / "omnibenchmark.yaml"
    yaml_file.write_text("""entrypoints:
  default: run.sh
""")

    entrypoint = _read_entrypoint(module_dir, "test_module")
    assert entrypoint == "run.sh"


@pytest.mark.short
def test_read_entrypoint_from_omnibenchmark_yaml_r_script(tmp_path):
    """Test reading R script entrypoint from omnibenchmark.yaml."""
    module_dir = tmp_path / "test_module"
    module_dir.mkdir()

    yaml_file = module_dir / "omnibenchmark.yaml"
    yaml_file.write_text("""entrypoints:
  default: run.R
""")

    entrypoint = _read_entrypoint(module_dir, "test_module")
    assert entrypoint == "run.R"


@pytest.mark.short
def test_read_entrypoint_from_config_cfg_with_warning(tmp_path):
    """Test reading entrypoint from config.cfg triggers FutureWarning."""
    module_dir = tmp_path / "test_module"
    module_dir.mkdir()

    # Create old-style config.cfg
    config_file = module_dir / "config.cfg"
    config_file.write_text("""[DEFAULT]
SCRIPT = run.py
""")

    # Should trigger FutureWarning
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        entrypoint = _read_entrypoint(module_dir, "test_module")

        assert entrypoint == "run.py"
        assert len(w) == 1
        assert issubclass(w[0].category, FutureWarning)
        assert "deprecated config.cfg" in str(w[0].message)
        assert "omnibenchmark.yaml" in str(w[0].message)


@pytest.mark.short
def test_read_entrypoint_prefers_omnibenchmark_yaml_over_config_cfg(tmp_path):
    """Test that omnibenchmark.yaml is preferred when both files exist."""
    module_dir = tmp_path / "test_module"
    module_dir.mkdir()

    # Create both files with different entrypoints
    yaml_file = module_dir / "omnibenchmark.yaml"
    yaml_file.write_text("""entrypoints:
  default: new_entrypoint.py
""")

    config_file = module_dir / "config.cfg"
    config_file.write_text("""[DEFAULT]
SCRIPT = old_entrypoint.py
""")

    # Should use omnibenchmark.yaml without warning
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        entrypoint = _read_entrypoint(module_dir, "test_module")

        assert entrypoint == "new_entrypoint.py"
        assert len(w) == 0  # No warning since we're using the new format


@pytest.mark.short
def test_read_entrypoint_no_config_file_returns_none(tmp_path):
    """Test that None is returned when no config file exists."""
    module_dir = tmp_path / "test_module"
    module_dir.mkdir()

    entrypoint = _read_entrypoint(module_dir, "test_module")
    assert entrypoint is None


@pytest.mark.short
def test_read_entrypoint_invalid_yaml_format_returns_none(tmp_path):
    """Test that invalid YAML format returns None."""
    module_dir = tmp_path / "test_module"
    module_dir.mkdir()

    yaml_file = module_dir / "omnibenchmark.yaml"
    yaml_file.write_text("""# Missing entrypoints key
something_else: value
""")

    entrypoint = _read_entrypoint(module_dir, "test_module")
    assert entrypoint is None


@pytest.mark.short
def test_read_entrypoint_yaml_missing_default_returns_none(tmp_path):
    """Test that YAML without 'default' key returns None."""
    module_dir = tmp_path / "test_module"
    module_dir.mkdir()

    yaml_file = module_dir / "omnibenchmark.yaml"
    yaml_file.write_text("""entrypoints:
  some_other_key: run.py
""")

    entrypoint = _read_entrypoint(module_dir, "test_module")
    assert entrypoint is None


@pytest.mark.short
def test_read_entrypoint_malformed_yaml_returns_none(tmp_path):
    """Test that malformed YAML returns None."""
    module_dir = tmp_path / "test_module"
    module_dir.mkdir()

    yaml_file = module_dir / "omnibenchmark.yaml"
    yaml_file.write_text("""entrypoints: [
  this is not valid yaml: {{{}
""")

    entrypoint = _read_entrypoint(module_dir, "test_module")
    assert entrypoint is None


@pytest.mark.short
def test_read_entrypoint_invalid_config_cfg_format_returns_none(tmp_path):
    """Test that invalid config.cfg format returns None."""
    module_dir = tmp_path / "test_module"
    module_dir.mkdir()

    config_file = module_dir / "config.cfg"
    config_file.write_text("""[DEFAULT]
# Missing SCRIPT key
OTHER_KEY = value
""")

    with warnings.catch_warnings(record=True):
        warnings.simplefilter("always")
        entrypoint = _read_entrypoint(module_dir, "test_module")
        assert entrypoint is None


@pytest.mark.short
def test_read_entrypoint_yaml_with_complex_structure(tmp_path):
    """Test reading entrypoint from YAML with other fields (future-proofing)."""
    module_dir = tmp_path / "test_module"
    module_dir.mkdir()

    # YAML with additional fields that might be added in the future
    yaml_file = module_dir / "omnibenchmark.yaml"
    yaml_file.write_text("""# OmniBenchmark Module Configuration

entrypoints:
  default: run.py

# Future fields (ignored for now)
module:
  name: my-module
  version: "1.0.0"
""")

    entrypoint = _read_entrypoint(module_dir, "test_module")
    assert entrypoint == "run.py"
