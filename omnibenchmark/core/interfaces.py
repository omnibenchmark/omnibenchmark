"""Protocol interfaces for artifact validation.

Protocols that decouple validation from model dependencies.
"""

from pathlib import Path
from typing import Protocol, Iterator


class ValidatableModule(Protocol):
    """Minimal interface for a module that can be validated."""

    @property
    def name(self) -> str:
        """Module identifier."""
        ...

    @property
    def path(self) -> Path:
        """Module directory containing validation.yaml."""
        ...

    @property
    def artifact_paths(self) -> dict[str, Path]:
        """Mapping of artifact names to file paths."""
        ...


class ValidatableBenchmark(Protocol):
    """Minimal interface for a benchmark that can be validated."""

    @property
    def path(self) -> Path:
        """Benchmark directory."""
        ...

    def __iter__(self) -> Iterator[ValidatableModule]:
        """Iterate over modules."""
        ...
