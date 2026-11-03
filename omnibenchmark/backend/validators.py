"""Run output validators: resolve the validators env and drive the generated
Snakefile (see ``omnibenchmark.core.validators`` for the pure parts)."""

import json
import shutil
import subprocess
import sys
from pathlib import Path

from omnibenchmark.core.validators import (
    VALIDATION_DIR,
    build_targets,
    declared_validators,
    verdict_status,
    write_validators_snakefile,
)
from omnibenchmark.logging import logger


def _resolve_env_directive(model, out_dir: Path, benchmark_dir: Path, env_id: str):
    """Return ``(directive_line, snakemake_backend_flag)`` for the validators env.

    Reuses the module resolver + the same directive mapping as the main run, so
    conda/apptainer/envmodules are activated identically. Host backend -> ("", None).
    """
    from omnibenchmark.backend.resolver import ModuleResolver
    from omnibenchmark.model.benchmark import SoftwareBackendEnum

    backend = model.get_software_backend()
    resolver = ModuleResolver(
        work_base_dir=out_dir / ".modules",
        output_dir=out_dir,
        software_backend=backend,
        software_environments=model.get_software_environments(),
        benchmark_dir=benchmark_dir,
    )
    resolved = resolver.resolve_environment(env_id, "validators")

    def _abs(ref: str) -> str:
        # conda yamls / local SIFs are resolved relative to out_dir, but the
        # validators Snakefile lives in out_dir/.validation — snakemake would
        # resolve a relative path against the wrong dir. Make file refs absolute.
        if "://" in ref or Path(ref).is_absolute():
            return ref
        return str((out_dir / ref).resolve())

    directive = ""
    if resolved is not None:
        bt = resolved.backend_type.value
        if bt == "conda":
            directive = f'conda: "{_abs(resolved.reference)}"'
        elif bt in ("apptainer", "docker"):
            directive = f'container: "{_abs(resolved.reference)}"'
        elif bt == "envmodules":
            directive = f'envmodules: "{resolved.reference}"'

    flag = {
        SoftwareBackendEnum.conda: "--use-conda",
        SoftwareBackendEnum.apptainer: "--use-singularity",
        SoftwareBackendEnum.docker: "--use-singularity",
        SoftwareBackendEnum.envmodules: "--use-envmodules",
    }.get(backend)
    return directive, flag


def run_validators(
    benchmark: str,
    out_dir: str = "out",
    cores: int = 1,
    force: bool = False,
) -> int:
    """Collect the declared validators, generate the Snakefile, run it.

    Returns non-zero if snakemake fails or any validator rejects its output.
    """
    from omnibenchmark.core.status.status import prepare_status
    from omnibenchmark.model import Benchmark

    benchmark_path = Path(benchmark)
    model = Benchmark.from_yaml(benchmark_path)
    cfg = model.get_validators_config()
    env = cfg.env if cfg else None
    if env and env not in model.get_software_environments():
        logger.error(f"validators env '{env}' is not defined in software_environments.")
        return 1

    out_path = Path(out_dir)
    if not out_path.is_dir():
        logger.error(f"Output dir '{out_path}' not found — run the benchmark first.")
        return 1

    benchmark_dir = benchmark_path.parent.absolute()
    validators = declared_validators(model, benchmark_dir)
    if not validators:
        logger.error("The plan declares no validators on any stage output.")
        return 1

    # Each node resolves its outputs in declaration order, so pairing them with
    # the output ids gives the file for every (stage, output id) without
    # matching filenames.
    _, _, exec_path_dict = prepare_status(
        str(benchmark_path),
        str(out_path),
        return_all=True,
        cache_dir=Path(".snakemake") / "repos",
    )
    produced = {
        (st.stage_id, output_id, f.output_file_path)
        for ep in exec_path_dict.values()
        for st in ep.exec_path.values()
        for output_id, f in zip(st.node.get_outputs(), st.output_stage_files)
        if f.output_file.exists()
    }
    targets = build_targets(produced, validators, out_path)
    if not targets:
        logger.error(
            "No produced outputs matched a validator; run the benchmark first."
        )
        return 1

    directive, backend_flag = (
        _resolve_env_directive(model, out_path, benchmark_dir, env)
        if env
        else ("", None)
    )
    snakefile = out_path / VALIDATION_DIR / "Snakefile"
    write_validators_snakefile(targets, directive, snakefile)

    bin_dir = Path(sys.executable).parent
    candidate = bin_dir / "snakemake"
    snakemake_bin = (
        str(candidate)
        if candidate.exists()
        else (shutil.which("snakemake") or "snakemake")
    )
    cmd = [
        snakemake_bin,
        "--snakefile",
        f"{VALIDATION_DIR}/Snakefile",
        "--cores",
        str(cores),
    ]
    if backend_flag:
        cmd.append(backend_flag)
    if force:
        cmd.append("--forceall")

    log_file = out_path / VALIDATION_DIR / "snakemake.log"
    logger.info(f"Validating {len(targets)} output file(s)")
    with open(log_file, "w") as f:
        result = subprocess.run(cmd, cwd=out_path, stdout=f, stderr=subprocess.STDOUT)
    if result.returncode != 0:
        tail = log_file.read_text(errors="replace").splitlines()[-20:]
        logger.error(
            f"Snakemake failed (exit {result.returncode}); full log: {log_file}\n"
            + "\n".join(tail)
        )
        return result.returncode

    # The rule records a failing validator instead of failing itself, so read
    # the results back to make the command usable as a CI gate.
    bad = 0
    for key, (rel, _, script, _) in targets.items():
        result_json = out_path / VALIDATION_DIR / f"{key}.json"
        try:
            status = verdict_status(json.loads(result_json.read_text())["exit_code"])
        except (OSError, ValueError, KeyError):
            status = "error"
        if status == "pass":
            continue
        bad += 1
        what = "Validation failed" if status == "fail" else "Validator could not run"
        logger.error(
            f"{what}: {rel} by {Path(script).name} "
            f"(log: {out_path / VALIDATION_DIR / key}.log)"
        )
    return 1 if bad else 0
