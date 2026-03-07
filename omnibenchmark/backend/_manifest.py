"""
Runtime manifest writer for omnibenchmark runs.

Writes `<output_dir>/.metadata/manifest.json` with a stable run UUID and
host metadata at the start of every `ob run` invocation.

See design/007-output-layout.md for the full schema specification and
limitations discussion.
"""

import json
import platform
import sys
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional


def write_run_manifest(
    output_dir: Path,
    run_id: Optional[str] = None,
) -> dict:
    """Write out/.metadata/manifest.json with a stable run UUID and host metadata.

    The *run_id* is the same UUID that the telemetry emitter uses as its
    OTLP trace_id, so every piece of provenance (telemetry.jsonl, manifest.json,
    log files) can be correlated by that single identifier.  When telemetry is
    disabled a fresh UUID is generated here.

    Known limitation: only the main executor node is documented here.  In a
    distributed Snakemake execution (cluster, cloud), worker nodes run on
    separate machines whose hardware is not captured.  See design
    007 §6 (Limitations) for a discussion and the proposed delta-update
    workaround.

    Fields written:
      run_id        – UUID4 string (from telemetry trace_id when available)
      ob_version    – omnibenchmark package version (from importlib.metadata)
      snakemake_cmd – exact snakemake argv list (patched in by _run_snakemake)
      timestamp     – ISO-8601 UTC timestamp of when the manifest was written
      hostname      – machine hostname
      platform      – OS platform string (e.g. "linux")
      os            – full OS release string (e.g. "Linux 6.x #1 SMP …")
      kernel        – kernel release (uname -r equivalent)
      cpu_count     – logical CPU count visible to the process
      cpu_model     – CPU model string where available (Linux /proc/cpuinfo)
      memory_total_mb  – total physical RAM in MiB (where available)
      python_version   – Python version string (e.g. "3.11.8")
      python_executable – path to the Python interpreter used to run the program
      gpu_devices      – list of NVIDIA GPU dicts {index, name, memory_total_mb}
                         from nvidia-smi; null if nvidia-smi is absent or fails
    """

    metadata_dir = output_dir / ".metadata"
    metadata_dir.mkdir(parents=True, exist_ok=True)

    if run_id is None:
        run_id = str(uuid.uuid4())

    uname = platform.uname()

    cpu_count = None
    try:
        import os

        cpu_count = os.cpu_count()
    except Exception:
        pass

    cpu_model = None
    try:
        if sys.platform.startswith("linux"):
            with open("/proc/cpuinfo") as fh:
                for line in fh:
                    if line.startswith("model name"):
                        cpu_model = line.split(":", 1)[1].strip()
                        break
        elif sys.platform == "darwin":
            import subprocess

            result = subprocess.run(
                ["sysctl", "-n", "machdep.cpu.brand_string"],
                capture_output=True,
                text=True,
                timeout=2,
            )
            if result.returncode == 0:
                cpu_model = result.stdout.strip() or None
    except Exception:
        pass

    memory_total_mb = None
    try:
        if sys.platform.startswith("linux"):
            with open("/proc/meminfo") as fh:
                for line in fh:
                    if line.startswith("MemTotal:"):
                        kb = int(line.split()[1])
                        memory_total_mb = kb // 1024
                        break
        elif sys.platform == "darwin":
            import subprocess

            result = subprocess.run(
                ["sysctl", "-n", "hw.memsize"],
                capture_output=True,
                text=True,
                timeout=2,
            )
            if result.returncode == 0:
                memory_total_mb = int(result.stdout.strip()) // (1024 * 1024)
    except Exception:
        pass

    # GPU devices via nvidia-smi (NVIDIA); optional.
    # nvidia-smi --query-gpu=... --format=csv,noheader,nounits is supported on
    # all modern driver versions and needs no extra Python dependencies.
    gpu_devices = None
    try:
        import subprocess as _sp

        result = _sp.run(
            [
                "nvidia-smi",
                "--query-gpu=index,name,memory.total",
                "--format=csv,noheader,nounits",
            ],
            capture_output=True,
            text=True,
            timeout=5,
        )
        if result.returncode == 0 and result.stdout.strip():
            gpu_devices = []
            for line in result.stdout.strip().splitlines():
                parts = [p.strip() for p in line.split(",")]
                entry = {"index": int(parts[0]), "name": parts[1]}
                if len(parts) >= 3:
                    try:
                        entry["memory_total_mb"] = int(parts[2])
                    except ValueError:
                        pass
                gpu_devices.append(entry)
    except Exception:
        pass

    try:
        from importlib.metadata import version as _pkg_version

        ob_version = _pkg_version("omnibenchmark")
    except Exception:
        ob_version = None

    manifest = {
        "run_id": run_id,
        "ob_version": ob_version,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "hostname": uname.node,
        "platform": sys.platform,
        "os": f"{uname.system} {uname.release} {uname.version}".strip(),
        "kernel": uname.release,
        "cpu_count": cpu_count,
        "cpu_model": cpu_model,
        "memory_total_mb": memory_total_mb,
        "python_version": platform.python_version(),
        "python_executable": sys.executable,
        "gpu_devices": gpu_devices,
    }

    manifest_path = metadata_dir / "manifest.json"
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)
        fh.write("\n")

    return manifest
