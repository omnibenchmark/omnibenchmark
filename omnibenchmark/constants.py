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
# `.metadata` is deliberately absent: it holds the run's provenance (benchmark
# YAML copy, modules.txt, manifest.json) and belongs in the archive.
INTERNAL_OUT_DIRS = frozenset({".modules", ".snakemake", ".envs", ".logs", ".cache"})
