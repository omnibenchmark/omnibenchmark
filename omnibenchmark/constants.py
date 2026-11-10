from enum import Enum


class LayoutDesign(Enum):
    Hierarchical = 1
    Spring = 2


# Archive compression
COMPRESSION_GZIP = "gzip"
DEFAULT_COMPRESSION = COMPRESSION_GZIP


# Directories under the results dir that hold execution state, not benchmark
# results: module checkouts, Snakemake bookkeeping (locks, logs, per-rule
# metadata), resolved environment files (which may symlink multi-GB images) and
# tool caches. They can reach many gigabytes, so nothing that walks `out/`
# looking for results should descend into them.
#
# `by-params` is the derived readable view of a flat tree (007 §3.1.2): symlinks
# into the results, so walking it would count every result twice.
#
# `.validation` holds `ob validate outputs` state: verdicts, logs and the
# generated Snakefile. The status report reads it in place.
#
# `.metadata` is deliberately absent: it holds the run's provenance (benchmark
# YAML copy, modules.txt, manifest.json) and belongs in the archive.
# `ob describe status --html` writes here by default; not a result.
STATUS_REPORT_FILE = "status_report.html"

INTERNAL_OUT_DIRS = frozenset(
    {".modules", ".snakemake", ".envs", ".logs", ".cache", ".validation", "by-params"}
)
