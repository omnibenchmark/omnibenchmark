"""Short unit tests for cli/progress.py pure functions and lightweight classes."""

import pytest

from omnibenchmark.cli.progress import (
    extract_rule_display_name,
    IncrementalLogReader,
    KeyboardListener,
)


@pytest.mark.short
class TestExtractRuleDisplayName:
    def test_all_rule(self):
        assert extract_rule_display_name("all") == "all"

    def test_default_rule(self):
        assert extract_rule_display_name("default") == "default"

    def test_empty_string(self):
        assert extract_rule_display_name("") == ""

    def test_three_part_name(self):
        # stageid_moduleid_paramhash → last two parts
        assert (
            extract_rule_display_name("methods_iris_method1_abc1234")
            == "method1_abc1234"
        )

    def test_two_part_name(self):
        assert extract_rule_display_name("stage_module") == "module"

    def test_one_part_name(self):
        assert extract_rule_display_name("singlepart") == "singlepart"

    def test_typical_rule(self):
        assert extract_rule_display_name("data_D1_d807c924") == "D1_d807c924"


@pytest.mark.short
class TestIncrementalLogReader:
    def test_read_new_lines_before_start(self, tmp_path):
        log = tmp_path / "test.log"
        log.write_text("line1\nline2\n")
        reader = IncrementalLogReader(log)
        # Not started — returns empty
        assert reader.read_new_lines() == []

    def test_start_creates_file_if_missing(self, tmp_path):
        log = tmp_path / "missing.log"
        reader = IncrementalLogReader(log)
        reader.start()
        assert log.exists()
        reader.stop()

    def test_reads_lines_after_start(self, tmp_path):
        log = tmp_path / "test.log"
        log.write_text("line1\nline2\n")
        reader = IncrementalLogReader(log)
        reader.start()
        lines = reader.read_new_lines()
        reader.stop()
        assert lines == ["line1", "line2"]

    def test_incremental_read(self, tmp_path):
        log = tmp_path / "test.log"
        log.write_text("first\n")
        reader = IncrementalLogReader(log)
        reader.start()
        reader.read_new_lines()
        # Append more content
        with open(log, "a") as f:
            f.write("second\n")
        lines = reader.read_new_lines()
        reader.stop()
        assert lines == ["second"]

    def test_get_tail_returns_last_n(self, tmp_path):
        log = tmp_path / "test.log"
        log.write_text("\n".join(f"line{i}" for i in range(10)) + "\n")
        reader = IncrementalLogReader(log, max_buffer_size=100)
        reader.start()
        reader.read_new_lines()
        tail = reader.get_tail(3)
        reader.stop()
        assert tail == ["line7", "line8", "line9"]

    def test_get_tail_fewer_than_n(self, tmp_path):
        log = tmp_path / "test.log"
        log.write_text("a\nb\n")
        reader = IncrementalLogReader(log)
        reader.start()
        reader.read_new_lines()
        tail = reader.get_tail(10)
        reader.stop()
        assert tail == ["a", "b"]

    def test_stop_idempotent(self, tmp_path):
        log = tmp_path / "test.log"
        reader = IncrementalLogReader(log)
        reader.stop()  # should not raise


@pytest.mark.short
class TestKeyboardListener:
    def test_get_key_when_disabled_returns_none(self):
        kb = KeyboardListener()
        # Not started → _enabled is False
        assert kb.get_key() is None

    def test_stop_when_not_started(self):
        kb = KeyboardListener()
        kb.stop()  # should not raise

    def test_start_in_non_tty_disables_listener(self):
        kb = KeyboardListener()
        kb.start()  # stdin is not a TTY in pytest → _enabled stays False
        assert not kb._enabled
        kb.stop()
