"""Run output validators: resolve the validators env and drive the generated
Snakefile (see ``omnibenchmark.core.validators`` for the pure parts)."""

import subprocess
from dataclasses import replace
from pathlib import Path

from omnibenchmark.backend.snakemake import (
    BACKEND_FLAGS,
    environment_directive,
    snakemake_bin,
)

from omnibenchmark.core.validators import (
    VALIDATION_DIR,
    build_targets,
    declared_validators,
    missing_scripts,
    read_verdict,
    stage_scripts,
    verdict_status,
    write_validators_snakefile,
)
from omnibenchmark.logging import logger


def _resolve_env_directive(model, out_dir: Path, benchmark_dir: Path, env_id: str):
    """``(directive, snakemake_flag)`` for the validators env, as the main run
    writes them. Host backend -> ``("", None)``; raises if the env cannot be
    resolved, so validators never fall back to the host silently."""
    from omnibenchmark.backend.resolver import ModuleResolver

    backend = model.get_software_backend()
    if backend not in BACKEND_FLAGS:
        return "", None
    resolver = ModuleResolver(
        work_base_dir=out_dir / ".modules",
        output_dir=out_dir,
        software_backend=backend,
        software_environments=model.get_software_environments(),
        benchmark_dir=benchmark_dir,
    )
    resolved = resolver.resolve_environment(env_id, "validators")
    if resolved is None:
        raise ValueError(f"could not resolve validators env '{env_id}'")
    # The main Snakefile sits in out/, the validators one in out/.validation/, and
    # snakemake reads env files relative to the Snakefile. Make the ones the
    # resolver placed under out/ absolute; image refs and module names stay.
    local = out_dir / resolved.reference
    if local.exists():
        resolved = replace(resolved, reference=str(local.resolve()))
    return environment_directive(resolved), BACKEND_FLAGS[backend]


def _drop_orphan_results(root: Path, targets) -> None:
    """Delete verdicts and logs of validators no longer in the plan."""
    keep = {root / "snakemake.log"}
    for key in targets:
        keep |= {root / f"{key}.json", root / f"{key}.log"}
    for p in root.rglob("*"):
        if p.suffix in (".json", ".log") and p.is_file() and p not in keep:
            p.unlink()


def run_validators(
    benchmark: str,
    out_dir: str = "out",
    cores: int = 1,
    force: bool = False,
) -> int:
    """Collect the declared validators, generate the Snakefile, run it.

    Returns non-zero if snakemake fails or any validator rejects its output.
    """
    from omnibenchmark.core.run_record import check_out_dir
    from omnibenchmark.core.status.status import prepare_status
    from omnibenchmark.model import Benchmark

    benchmark_path = Path(benchmark)
    model = Benchmark.from_yaml(benchmark_path)
    cfg = model.get_validators_config()
    # An undeclared env is rejected when the plan loads.
    env = cfg.env if cfg else None

    out_path = Path(out_dir)
    if not out_path.is_dir():
        logger.error(f"Output dir '{out_path}' not found — run the benchmark first.")
        return 1
    errors, notes = check_out_dir(model, out_path)
    for msg in errors:
        logger.error(msg)
    if errors:
        return 1
    for msg in notes:
        logger.warning(msg)

    benchmark_dir = benchmark_path.parent.absolute()
    validators = declared_validators(model, benchmark_dir)
    if not validators:
        logger.error("The plan declares no validators on any stage output.")
        return 1
    missing = missing_scripts(validators)
    if missing:
        logger.error(f"Validator script(s) not found: {', '.join(missing)}")
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

    try:
        directive, backend_flag = (
            _resolve_env_directive(model, out_path, benchmark_dir, env)
            if env
            else ("", None)
        )
    except ValueError as e:
        logger.error(str(e))
        return 1
    snakefile = out_path / VALIDATION_DIR / "Snakefile"
    write_validators_snakefile(targets, directive, snakefile)
    stage_scripts(targets, out_path)
    _drop_orphan_results(out_path / VALIDATION_DIR, targets)

    cmd = [
        snakemake_bin(),
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
        log = log_file.read_text(errors="replace")
        msg = f"Snakemake failed (exit {result.returncode}); full log: {log_file}\n"
        msg += "\n".join(log.splitlines()[-20:])
        if "cannot be locked" in log:
            msg += (
                f"\nAnother Snakemake run holds the lock on {out_path}, most "
                "likely `ob run`: wait for it to finish. If none is running, "
                f"run `snakemake --unlock` in {out_path}."
            )
        logger.error(msg)
        return result.returncode

    # The rule records a failing validator instead of failing itself, so read
    # the results back to make the command usable as a CI gate.
    bad = 0
    for key, (rel, _, script, _) in targets.items():
        verdict = read_verdict(out_path / VALIDATION_DIR / f"{key}.json")
        status = verdict_status(verdict.get("exit_code"))
        if status == "pass":
            continue
        bad += 1
        what = "Validation failed" if status == "fail" else "Validator could not run"
        logger.error(
            f"{what}: {rel} by {Path(script).name} "
            f"(log: {out_path / VALIDATION_DIR / key}.log)"
        )
    return 1 if bad else 0
