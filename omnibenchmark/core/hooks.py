"""Git-style run hooks (design 016): pre-run and post-run.

A hook is the command `hooks.<event>` in omnibenchmark.cfg, or the executable
`<hooks.path>/<event>` (`OB_HOOKS_PATH` overrides `hooks.path`). Nothing runs
unless one is configured. Hooks observe a run; they must not change it.
"""

import hashlib
import os
import shlex
import subprocess
from pathlib import Path

from omnibenchmark import config as _config
from omnibenchmark.logging import logger

EVENTS = ("pre-run", "post-run")
# Not results: the hooks' own files, built envs, module checkouts (GBs, 10^5 files).
SKIP = {".hooks", ".snakemake", ".modules"}
VETO = 3  # ob's exit code when pre-run stops the run


def find(event):
    """(argv, source) of the hook for event, or None."""
    cfg = _config.config
    cmd = cfg.get("hooks", event)
    if cmd:
        return shlex.split(cmd), f"hooks.{event}: {cmd}"
    d = os.environ.get("OB_HOOKS_PATH") or cfg.get("hooks", "path")
    if d:
        f = Path(d).expanduser() / event
        if f.is_file() and os.access(f, os.X_OK):
            return [str(f)], str(f)
    return None


def log_active():
    active = [f"{e} ({h[1]})" for e in EVENTS if (h := find(e))]
    if active:
        logger.info("Running with hooks: " + ", ".join(active))


def _snapshot(out_dir):
    """size and mtime of every file in out_dir outside SKIP."""
    snap = {}
    for root, dirs, files in os.walk(out_dir):
        if root == str(out_dir):
            dirs[:] = [d for d in dirs if d not in SKIP]
        for f in files:
            p = Path(root, f)
            try:
                s = p.stat()
            except OSError:  # gone, or a broken symlink
                continue
            snap[p] = (s.st_size, s.st_mtime_ns)
    return snap


def run(event, plan, out_dir, *args, env=None):
    """Run event's hook with <plan> <out_dir> *args; return its exit status (0 if none)."""
    hook = find(event)
    if hook is None:
        return 0
    plan, out_dir = Path(plan).resolve(), Path(out_dir).resolve()
    hooks_out = out_dir / ".hooks"
    if event == "post-run":  # pre-run writes nothing; the hook may mkdir it
        hooks_out.mkdir(parents=True, exist_ok=True)
    full_env = {
        **os.environ,
        **(env or {}),
        "OB_HOOK": event,
        "OB_HOOKS_OUT": str(hooks_out),
    }
    timeout = float(_config.config.get("hooks", "timeout", 60))

    before_plan = hashlib.sha256(plan.read_bytes()).digest()
    before_out = _snapshot(out_dir) if event == "post-run" else None
    try:
        rc = subprocess.run(
            hook[0] + [str(plan), str(out_dir), *map(str, args)],
            env=full_env,
            stdout=2,  # stderr's fd: sys.stderr may have none (embedded, CliRunner)
            timeout=timeout,
        ).returncode
    except subprocess.TimeoutExpired:
        logger.warning(f"{event} hook timed out after {timeout:g}s")
        rc = 1
    except OSError as e:
        logger.warning(f"{event} hook failed to start: {e}")
        rc = 1

    if hashlib.sha256(plan.read_bytes()).digest() != before_plan:
        logger.error(f"plan modified by {event} hook: {plan}")
        rc = rc or 1
    if before_out is not None:
        diff = _snapshot(out_dir).items() ^ before_out.items()
        if diff:
            changed = sorted({str(p.relative_to(out_dir)) for p, _ in diff})
            logger.warning(f"modified by post-run hook: {', '.join(changed)}")
    return rc
