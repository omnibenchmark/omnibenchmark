"""Modern progress display using rich."""

from rich.console import Console
from rich.progress import (
    Progress,
    SpinnerColumn,
    TextColumn,
    BarColumn,
    TaskProgressColumn,
    TimeElapsedColumn,
)
from rich.table import Table


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
            TaskProgressColumn(),
            TimeElapsedColumn(),
            console=self.console,
            refresh_per_second=10,  # Update 10 times per second for responsiveness
            transient=True,  # Progress bar disappears after completion
        )
        self._progress.start()
        self._task = self._progress.add_task(description, total=total)

    def update(self, advance: int = 1):
        """Update progress bar."""
        if self._progress and self._task is not None:
            self._progress.update(self._task, advance=advance)
            self._progress.refresh()  # Force immediate refresh

    def finish(self):
        """Finish and hide progress bar."""
        if self._progress:
            self._progress.stop()
            self._progress = None
            self._task = None

    def status(self, message: str):
        """Print a status message."""
        self.console.print(f"[dim]→[/dim] {message}")

    def success(self, message: str):
        """Print a success message."""
        self.console.print(f"[green]{message}[/green]")

    def error(self, message: str):
        """Print an error message."""
        self.console.print(f"[red]{message}[/red]")

    def section(self, title: str):
        """Print a section header."""
        self.console.print()
        self.console.rule(f"[bold cyan]{title}[/bold cyan]")

    def summary(self, data: dict):
        """Print a summary table."""
        table = Table(show_header=False, box=None, padding=(0, 2))
        table.add_column(style="dim")
        table.add_column()

        for key, value in data.items():
            table.add_row(f"{key}:", str(value))

        self.console.print(table)
