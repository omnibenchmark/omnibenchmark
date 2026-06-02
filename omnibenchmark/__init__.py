"""
omnibenchmark — define, run, and share reproducible benchmarks.

Package layering (imports flow top → bottom; nothing imports upward):

  cli            user-facing ``ob`` commands
  ────────────────────────────────────────────────────────
  core           execution engine: build the DAG, run a benchmark
  archive  remote  versioning  backend   higher-level services
  git            two-tier repository cache
  ────────────────────────────────────────────────────────
  model          benchmark definition (Pydantic) + validation
  dag            generic DAG primitives
  config  constants  logging  progress  error_formatting   leaf utilities

The only remaining import cycle is ``archive <-> remote``; everything else
forms a clean DAG.
"""

try:
    from importlib.metadata import version

    __version__ = version("omnibenchmark")

except Exception:
    __version__ = "0.0.0"  # Fallback version for development
