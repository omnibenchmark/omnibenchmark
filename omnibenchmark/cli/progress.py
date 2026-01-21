"""Modern progress display using rich with interactive log viewing."""

import sys
import select
import threading
import time
from pathlib import Path
from typing import Optional
from collections import deque

from rich.console import Console, RenderableType
from rich.live import Live
from rich.panel import Panel
from rich.progress import (
    Progress,
    SpinnerColumn,
    TextColumn,
    BarColumn,
    MofNCompleteColumn,
    TimeElapsedColumn,
)
from rich.table import Table
from rich.text import Text


def extract_rule_display_name(rule_name: str) -> str:
    """
    Extract a clean display name from a snakemake rule name.

    Rule names follow pattern: stageid_moduleid_paramhash
    e.g., "methods_iris_method1_abc1234" -> "method1_abc1234"
         "metrics_iris_report_def5678" -> "report_def5678"
         "all" -> "all"
    """
    if not rule_name or rule_name in ("all", "default"):
        return rule_name

    parts = rule_name.split("_")
    if len(parts) >= 3:
        return "_".join(parts[-2:])
    elif len(parts) == 2:
        return parts[-1]
    return rule_name


class ProgressDisplay:
    """Clean progress display with rich."""

    def __init__(self):
        self.console = Console()
        self._progress = None
        self._task = None

    def start_task(self, description: str, total: int):
        """Start a progress bar for a task."""
        self._progress = Progress(
            SpinnerColumn(),
            TextColumn("[bold blue]{task.description}"),
            BarColumn(),
            MofNCompleteColumn(),
            TimeElapsedColumn(),
            console=self.console,
            refresh_per_second=10,
            transient=True,
        )
        self._progress.start()
        self._task = self._progress.add_task(description, total=total)

    def update(self, advance: int = 0, description: str = None):
        """Update progress bar."""
        if self._progress and self._task is not None:
            if description:
                self._progress.update(
                    self._task, advance=advance, description=description
                )
            else:
                self._progress.update(self._task, advance=advance)
            self._progress.refresh()

    def finish(self):
        """Finish and hide progress bar."""
        if self._progress:
            self._progress.stop()
            self._progress = None
            self._task = None

    def status(self, message: str):
        self.console.print(f"[dim]→[/dim] {message}")

    def success(self, message: str):
        self.console.print(f"[green]{message}[/green]")

    def error(self, message: str):
        self.console.print(f"[red]{message}[/red]")

    def section(self, title: str):
        self.console.print()
        self.console.rule(f"[bold cyan]{title}[/bold cyan]")

    def summary(self, data: dict):
        table = Table(show_header=False, box=None, padding=(0, 2))
        table.add_column(style="dim")
        table.add_column()
        for key, value in data.items():
            table.add_row(f"{key}:", str(value))
        self.console.print(table)


class KeyboardListener:
    """Non-blocking keyboard listener."""

    def __init__(self):
        self._old_settings = None
        self._fd = None
        self._enabled = False

    def start(self):
        try:
            import termios
            import tty

            if not sys.stdin.isatty():
                return
            self._fd = sys.stdin.fileno()
            self._old_settings = termios.tcgetattr(self._fd)
            tty.setcbreak(self._fd)
            self._enabled = True
        except (ImportError, termios.error, AttributeError, OSError):
            pass

    def stop(self):
        if self._old_settings is not None and self._fd is not None:
            try:
                import termios

                termios.tcsetattr(self._fd, termios.TCSADRAIN, self._old_settings)
            except (ImportError, termios.error, AttributeError, OSError):
                pass
        self._enabled = False
        self._old_settings = None
        self._fd = None

    def get_key(self) -> Optional[str]:
        if not self._enabled:
            return None
        try:
            if select.select([sys.stdin], [], [], 0)[0]:
                return sys.stdin.read(1)
        except (OSError, ValueError):
            pass
        return None


class IncrementalLogReader:
    """Efficiently reads new log lines as they're appended to a file."""

    def __init__(self, log_file: Path, max_buffer_size: int = 100):
        self.log_file = Path(log_file)
        self.max_buffer_size = max_buffer_size
        self._file_handle = None
        self._buffer = deque(maxlen=max_buffer_size)
        self._position = 0

    def start(self):
        """Open the log file for reading."""
        try:
            if not self.log_file.exists():
                # Create empty file if it doesn't exist
                self.log_file.touch()
            self._file_handle = open(self.log_file, "r", errors="replace")
            self._position = 0
        except Exception:
            self._file_handle = None

    def read_new_lines(self) -> list[str]:
        """Read any new lines added since last read."""
        if self._file_handle is None:
            return []

        try:
            # Seek to last known position
            self._file_handle.seek(self._position)
            new_lines = self._file_handle.readlines()

            # Update position
            self._position = self._file_handle.tell()

            # Add to buffer
            for line in new_lines:
                self._buffer.append(line.rstrip())

            return [line.rstrip() for line in new_lines]
        except Exception:
            return []

    def get_tail(self, n: int = 20) -> list[str]:
        """Get the last n lines from the buffer."""
        if n >= len(self._buffer):
            return list(self._buffer)
        return list(self._buffer)[-n:]

    def stop(self):
        """Close the file handle."""
        if self._file_handle:
            try:
                self._file_handle.close()
            except Exception:
                pass
            self._file_handle = None


class InteractiveProgress:
    """
    Interactive progress with keyboard toggle for live log viewing.

    Shows: [spinner] description [bar] 3/10 rule_name 0:01:23 [f=logs]
    Press 'f' to toggle between progress and log view.
    """

    def __init__(self, log_file: Path, tail_lines: int = 30):
        self.console = Console()
        self.log_file = Path(log_file)
        self.tail_lines = tail_lines

        self._show_log = False
        self._description = "Starting..."
        self._completed = 0
        self._total = 0
        self._failed_rules: list[str] = []
        self._current_rule: Optional[str] = None

        # Scroll state for log view
        self._scroll_offset = 0  # 0 = bottom (latest), positive = scrolled up
        self._auto_scroll = True  # Auto-scroll to bottom when new lines arrive

        # Progress bar - always created, just shown/hidden
        self._progress: Optional[Progress] = None
        self._task = None

        # Live display - manages rendering of both views
        self._live: Optional[Live] = None

        # Efficient log reader with buffer
        self._log_reader = IncrementalLogReader(log_file, max_buffer_size=200)

        self._keyboard = KeyboardListener()

        # Background keyboard monitoring
        self._keyboard_thread: Optional[threading.Thread] = None
        self._stop_keyboard_thread = threading.Event()
        self._view_changed = threading.Event()

    def start(self, description: str, total: int):
        self._description = description
        self._total = total
        self._completed = 0

        self._keyboard.start()
        self._log_reader.start()

        # Create progress bar (but don't start it yet)
        self._progress = Progress(
            SpinnerColumn(),
            TextColumn("[bold blue]{task.description}"),
            BarColumn(bar_width=30),
            MofNCompleteColumn(),
            TextColumn("{task.fields[rule]}"),
            TimeElapsedColumn(),
            TextColumn("[dim]f=logs, q=back[/dim]"),
            console=self.console,
            refresh_per_second=10,
            transient=True,
        )

        rule_text = f"[dim]{self._current_rule}[/dim]" if self._current_rule else ""
        self._task = self._progress.add_task(
            self._description,
            total=self._total,
            completed=self._completed,
            rule=rule_text,
        )

        # Start unified live display with higher refresh rate
        self._live = Live(
            self._render_current_view(),
            console=self.console,
            refresh_per_second=15,  # Higher refresh for snappier feel
            transient=True,
        )
        self._live.start()

        # Start background keyboard monitoring thread
        self._stop_keyboard_thread.clear()
        self._keyboard_thread = threading.Thread(
            target=self._keyboard_monitor_loop, daemon=True
        )
        self._keyboard_thread.start()

    def _render_current_view(self) -> RenderableType:
        """Render either progress or log view based on current mode."""
        if self._show_log:
            return self._render_log_view()
        else:
            return self._render_progress_view()

    def _render_progress_view(self) -> RenderableType:
        """Render the progress bar view."""
        # Update progress bar state
        if self._task is not None:
            rule_text = f"[dim]{self._current_rule}[/dim]" if self._current_rule else ""
            self._progress.update(
                self._task,
                completed=self._completed,
                description=self._description,
                rule=rule_text,
            )
        return self._progress

    def _render_log_view(self) -> RenderableType:
        """Render the live log view with syntax highlighting and scrolling."""
        # Read any new lines (efficient - only reads what's new)
        new_lines = self._log_reader.read_new_lines()

        # Auto-scroll to bottom if enabled and new content arrived
        if self._auto_scroll and new_lines:
            self._scroll_offset = 0

        # Get all available lines from buffer
        all_lines = list(self._log_reader._buffer)
        total_lines = len(all_lines)

        # Calculate which lines to show based on scroll offset
        # offset=0 means show the latest lines (bottom)
        # offset>0 means scroll up by that many lines
        if self._scroll_offset == 0:
            # Show the last N lines (tail)
            lines = (
                all_lines[-self.tail_lines :]
                if total_lines > self.tail_lines
                else all_lines
            )
            start_line_num = max(0, total_lines - len(lines))
        else:
            # Scrolled up - show lines from offset
            end_idx = total_lines - self._scroll_offset
            start_idx = max(0, end_idx - self.tail_lines)
            lines = all_lines[start_idx:end_idx]
            start_line_num = start_idx

        # Clamp scroll offset to valid range
        max_offset = max(0, total_lines - self.tail_lines)
        if self._scroll_offset > max_offset:
            self._scroll_offset = max_offset

        formatted = []
        for line in lines:
            if line.startswith("rule ") or line.startswith("localrule "):
                formatted.append(f"[cyan]{line}[/cyan]")
            elif "Error" in line or "error:" in line.lower():
                formatted.append(f"[red]{line}[/red]")
            elif "Finished job" in line:
                formatted.append(f"[green]{line}[/green]")
            elif line.strip().startswith("input:") or line.strip().startswith(
                "output:"
            ):
                formatted.append(f"[dim]{line}[/dim]")
            else:
                formatted.append(line)

        # Header with progress info and scroll indicators
        header = f"[bold]{self._completed}/{self._total}[/bold]"
        if self._current_rule:
            header += f"  [cyan]{self._current_rule}[/cyan]"
        if self._failed_rules:
            header += f"  [red]{len(self._failed_rules)} failed[/red]"

        # Scroll position indicator
        if total_lines > self.tail_lines:
            if self._scroll_offset > 0:
                header += f"  [yellow]↑ {self._scroll_offset} lines up[/yellow]"
            else:
                header += "  [dim]⬇ at bottom[/dim]"

        log_text = (
            "\n".join(formatted) if formatted else "[dim]Waiting for output...[/dim]"
        )

        # Footer with scroll indicators
        scroll_hints = []
        can_scroll_up = start_line_num > 0
        can_scroll_down = self._scroll_offset > 0

        if can_scroll_up:
            scroll_hints.append("↑/PgUp=scroll up")
        if can_scroll_down:
            scroll_hints.append("↓/PgDn=scroll down")
        if can_scroll_down:
            scroll_hints.append("End=bottom")
        if can_scroll_up:
            scroll_hints.append("Home=top")

        content = f"{header}\n[dim]{'─' * 50}[/dim]\n{log_text}"

        subtitle_parts = ["q/esc/p=back"]
        if scroll_hints:
            subtitle_parts.append(" | ".join(scroll_hints))
        subtitle = "[dim]" + " | ".join(subtitle_parts) + "[/dim]"

        return Panel(
            Text.from_markup(content),
            title="[bold]Live Log[/bold]",
            subtitle=subtitle,
            border_style="yellow",
        )

    def _keyboard_monitor_loop(self):
        """Background thread that continuously monitors keyboard input."""
        last_log_update = 0
        log_update_interval = (
            0.1  # Update log view every 100ms when in auto-scroll mode
        )

        while not self._stop_keyboard_thread.is_set():
            key = self._keyboard.get_key()
            if key:
                if key.lower() == "f" and not self._show_log:
                    # Switch to log view - just toggle flag
                    self._show_log = True
                    self._scroll_offset = 0
                    self._auto_scroll = True
                    self._view_changed.set()
                    # Force immediate update
                    if self._live:
                        self._live.update(self._render_current_view(), refresh=True)
                elif key.lower() in ("q", "p"):
                    if self._show_log:
                        # Switch back to progress view - instant toggle
                        self._show_log = False
                        self._view_changed.set()
                        # Force immediate update
                        if self._live:
                            self._live.update(self._render_current_view(), refresh=True)
                elif key == "\x1b":
                    # ESC sequence - might be arrow key or standalone ESC
                    # Try to read the rest of the escape sequence
                    next_key = self._keyboard.get_key()
                    if next_key == "[" and self._show_log:
                        # Arrow key sequence in log view
                        arrow_key = self._keyboard.get_key()
                        updated = False
                        if arrow_key == "A":  # Up arrow
                            self._scroll_offset = min(
                                self._scroll_offset + 1, len(self._log_reader._buffer)
                            )
                            self._auto_scroll = False
                            updated = True
                        elif arrow_key == "B":  # Down arrow
                            self._scroll_offset = max(0, self._scroll_offset - 1)
                            if self._scroll_offset == 0:
                                self._auto_scroll = True
                            updated = True
                        elif arrow_key == "5":  # PgUp
                            self._keyboard.get_key()  # consume '~'
                            self._scroll_offset = min(
                                self._scroll_offset + self.tail_lines,
                                len(self._log_reader._buffer),
                            )
                            self._auto_scroll = False
                            updated = True
                        elif arrow_key == "6":  # PgDn
                            self._keyboard.get_key()  # consume '~'
                            self._scroll_offset = max(
                                0, self._scroll_offset - self.tail_lines
                            )
                            if self._scroll_offset == 0:
                                self._auto_scroll = True
                            updated = True
                        elif arrow_key == "H":  # Home
                            self._scroll_offset = len(self._log_reader._buffer)
                            self._auto_scroll = False
                            updated = True
                        elif arrow_key == "F":  # End
                            self._scroll_offset = 0
                            self._auto_scroll = True
                            updated = True

                        if updated and self._live:
                            self._live.update(self._render_current_view(), refresh=True)
                    elif next_key is None and self._show_log:
                        # Standalone ESC - close log view
                        self._show_log = False
                        self._view_changed.set()
                        if self._live:
                            self._live.update(self._render_current_view(), refresh=True)

            # When in log view with auto-scroll, periodically update to show new lines
            current_time = time.time()
            if (
                self._show_log
                and self._auto_scroll
                and (current_time - last_log_update) >= log_update_interval
            ):
                last_log_update = current_time
                if self._live:
                    self._live.update(self._render_current_view(), refresh=True)

            # Small sleep to avoid busy-waiting (50ms = 20 checks/sec)
            time.sleep(0.05)

    def check_keyboard(self):
        """Check if view changed (for backwards compatibility)."""
        # This is now handled by background thread, but we keep the method
        # for any existing calls from run.py
        if self._view_changed.is_set():
            self._view_changed.clear()
            return True
        return False

    def update(
        self, advance: int = 0, description: str = None, current_rule: str = None
    ):
        """Update progress state and refresh display."""
        if advance:
            self._completed += advance
        if description:
            self._description = description
        if current_rule is not None:
            self._current_rule = extract_rule_display_name(current_rule)

        # Update the unified live display (automatically shows correct view)
        if self._live:
            self._live.update(self._render_current_view())

    def add_failed_rule(self, rule: str):
        self._failed_rules.append(rule)

    def finish(self):
        """Clean up all resources."""
        # Stop keyboard monitoring thread
        self._stop_keyboard_thread.set()
        if self._keyboard_thread and self._keyboard_thread.is_alive():
            self._keyboard_thread.join(timeout=0.5)

        self._keyboard.stop()
        self._log_reader.stop()
        if self._live:
            self._live.stop()
            self._live = None
        # Progress is managed by Live, no need to stop separately
        self._progress = None
        self._task = None

    @property
    def failed_rules(self) -> list[str]:
        return self._failed_rules.copy()

    @property
    def completed(self) -> int:
        return self._completed
