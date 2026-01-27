"""Tests for argument parsing utilities."""

from click.testing import CliRunner

from omnibenchmark.cli.run import run
from omnibenchmark.cli.utils.args import parse_extra_args


class TestParseExtraArgs:
    """Test the parse_extra_args function.

    Note: In Click, when using ignore_unknown_options=True and allow_extra_args=True,
    the -- separator is automatically stripped by Click, and ctx.args only contains
    arguments that were specified after --.
    """

    def test_empty_args_returns_empty_dict(self):
        """When ctx.args is empty (no args after --), should return empty dict."""
        args = []
        result = parse_extra_args(args)
        assert result == {}

    def test_flag_only(self):
        """Parse flag without value (Click already stripped --)."""
        args = ["--dryrun"]
        result = parse_extra_args(args)
        assert result == {"dryrun": True}

    def test_flag_and_value(self):
        """Parse flag with single value (Click already stripped --)."""
        args = ["--cores", "4"]
        result = parse_extra_args(args)
        assert result == {"cores": "4"}

    def test_multiple_flags(self):
        """Parse multiple flags (Click already stripped --)."""
        args = ["--cores", "4", "--dryrun", "--verbose"]
        result = parse_extra_args(args)
        assert result == {"cores": "4", "dryrun": True, "verbose": True}

    def test_flag_and_multiple_values(self):
        """Parse flag with multiple values (Click already stripped --)."""
        args = ["--set-threads", "rule1=4", "rule2=8"]
        result = parse_extra_args(args)
        assert result == {"set-threads": ["rule1=4", "rule2=8"]}

    def test_complex_example(self):
        """Test a complex real-world example (Click already stripped --)."""
        args = [
            "--dryrun",
            "--verbose",
            "--set-threads",
            "rule1=4",
            "rule2=8",
            "--configfile",
            "config.yaml",
        ]
        result = parse_extra_args(args)
        assert result == {
            "dryrun": True,
            "verbose": True,
            "set-threads": ["rule1=4", "rule2=8"],
            "configfile": "config.yaml",
        }

    def test_mixed_flags_and_values(self):
        """Test various combinations of flags and values."""
        args = [
            "--jobs",
            "6",
            "--default-resources",
            "mem_mb=4000",
            "runtime=600",
            "--printshellcmds",
        ]
        result = parse_extra_args(args)
        assert result == {
            "jobs": "6",
            "default-resources": ["mem_mb=4000", "runtime=600"],
            "printshellcmds": True,
        }


class TestCLIUnknownOptionsError:
    """Test that unknown options without -- separator produce clear errors."""

    def test_unknown_option_without_separator_errors(self):
        """Unknown options without -- should error, not silently pass through."""
        runner = CliRunner()
        result = runner.invoke(
            run,
            ["bench.yaml", "--dryrun"],  # --dryrun is unknown to ob
        )
        assert result.exit_code != 0
        assert "No such option: --dryrun" in result.output

    def test_unknown_option_with_separator_succeeds(self, tmp_path):
        """Unknown options after -- should be accepted."""
        runner = CliRunner()
        # Just test that --dryrun after -- doesn't cause "No such option" error
        result = runner.invoke(
            run,
            ["nonexistent.yaml", "--", "--dryrun"],
        )
        # Will fail for other reasons (file not found), but NOT "No such option"
        assert "No such option: --dryrun" not in result.output
