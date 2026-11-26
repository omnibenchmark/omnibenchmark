"""A node that did not produce its outputs (design 014).

The node wrapper (`backend/ob_node.sh`) writes `failure.json` into the node
directory; everything in Python reads it through `Failure`.
"""

import json
import os
import re
from collections import defaultdict
from dataclasses import asdict, dataclass, field
from enum import Enum
from pathlib import Path
from typing import Iterable, List, Optional, Tuple, Union

FAILURE_FILE = "failure.json"


class FailureKind(str, Enum):
    FAILED_TO_START = "failed_to_start"  # entrypoint, interpreter or package missing
    CRASHED = "crashed"  # ran and exited non-zero, or was killed
    TIMED_OUT = "timed_out"  # killed by ob at the timeout
    MISSING_OUTPUT = "missing_output"  # exit 0, a declared output not written
    UPSTREAM = "upstream"  # not run: an input comes from a failed node


@dataclass(frozen=True)
class Failure:
    kind: FailureKind
    reason: str
    exit_code: Optional[int] = None  # None when the module was not run
    elapsed_s: Optional[int] = None
    timeout_s: Optional[int] = None
    upstream: Tuple[str, ...] = ()
    log: Optional[str] = None
    node_dir: str = field(default="", compare=False)  # from the loader

    @property
    def signal(self) -> Optional[int]:
        """The signal that killed the module (exit codes above 128)."""
        return self.exit_code - 128 if self.exit_code and self.exit_code > 128 else None

    @classmethod
    def load(
        cls, node_dir: Union[str, Path], out_dir: Path = Path(".")
    ) -> Optional["Failure"]:
        """The node's failure, or None if it has no record."""
        path = Path(out_dir) / node_dir / FAILURE_FILE
        if not path.exists():
            return None
        record = json.loads(path.read_text())
        record["kind"] = FailureKind(record["kind"])
        record["upstream"] = tuple(record.get("upstream") or ())
        return cls(**record, node_dir=str(node_dir))

    def to_dict(self) -> dict:
        record = asdict(self)
        record.pop("node_dir")
        record["kind"] = self.kind.value
        record["upstream"] = list(self.upstream)
        return record


def parse_duration(text: str) -> int:
    """`90s`, `30m`, `2h`, `1d` to seconds. The unit is required: `runtime`
    in resources is in minutes, so a bare number would be read both ways."""
    m = re.fullmatch(r"\s*(\d+)\s*([smhd])\s*", text)
    if not m or int(m[1]) == 0:
        raise ValueError(f"invalid duration {text!r}: expected e.g. 90s, 30m, 2h, 1d")
    return int(m[1]) * {"s": 1, "m": 60, "h": 3600, "d": 86400}[m[2]]


def node_dir(node) -> str:
    """Where a node's failure.json goes: its directory, else its first
    output's (the same fallback the gather manifest uses)."""
    return node.node_dir or os.path.dirname(next(iter(node.outputs.values()), ""))


def collect(out_dir: Path, nodes: Iterable) -> List[Failure]:
    """Every recorded failure among *nodes* (ResolvedNode), in plan order."""
    found = (Failure.load(node_dir(n), out_dir) for n in nodes)
    return [f for f in found if f is not None]


def clear(out_dir: Path, nodes: Iterable, kinds: Optional[set] = None) -> int:
    """Delete the record and outputs of failed nodes (`ob run --retry-failed`),
    so Snakemake runs them again. Returns how many were cleared."""
    cleared = 0
    for node in nodes:
        failure = Failure.load(node_dir(node), out_dir)
        if failure is None or (kinds and failure.kind.value not in kinds):
            continue
        for path in node.outputs.values():
            (Path(out_dir) / path).unlink(missing_ok=True)
        (Path(out_dir) / node_dir(node) / FAILURE_FILE).unlink()
        cleared += 1
    return cleared


def summary(failures: List[Failure], total: int) -> str:
    """One line per failure, grouped by kind."""
    by_kind = defaultdict(list)
    for f in failures:
        by_kind[f.kind].append(f)
    lines = [f"{len(failures)} of {total} nodes failed"]
    for kind in FailureKind:
        if by_kind[kind]:
            lines.append(f"  {kind.value} ({len(by_kind[kind])})")
            lines += [f"    {f.node_dir}   {f.reason}" for f in by_kind[kind]]
    return "\n".join(lines)
