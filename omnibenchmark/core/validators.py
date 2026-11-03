"""Output validators: discover validator scripts, generate the Snakefile that runs
them against the exact produced files, and read back per-output pass/fail results.

Validators are declared on a stage output (``validators: [{script}]``). Each
receives a single produced output file path and signals pass/fail by EXIT CODE
(0 = pass, 126/127 = could not run, anything else = fail); stdout/stderr are
captured to a per-output log.

Targets are the concrete files the benchmark produced (from the status machinery),
keyed by the output id they were produced for, so no filename matching or
filesystem globbing is needed.

Pure helpers; ``omnibenchmark.backend.validators`` runs them. The generated
Snakefile is run from the benchmark out/ dir; results and logs land under
``out/.validation/`` mirroring each output's path.
"""

import hashlib
import json
import os
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

# Interpreter by validator extension. Anything else falls back to bash.
INTERP = {".py": "python3", ".R": "Rscript", ".sh": "bash"}

VALIDATION_DIR = ".validation"

# Key: (stage_id, output_id). Value: [(interpreter, absolute script path), ...].
ValidatorMap = Dict[Tuple[str, str], List[Tuple[str, str]]]


def declared_validators(model, benchmark_dir: Path) -> ValidatorMap:
    """Map ``(stage_id, output_id)`` to the validators declared on that output."""
    return {
        (st.id, o.id): [
            (
                INTERP.get(Path(v.script).suffix, "bash"),
                str((benchmark_dir / v.script).resolve()),
            )
            for v in o.validators
        ]
        for st in model.stages
        for o in st.outputs
        if o.validators
    }


def missing_scripts(validators: ValidatorMap) -> List[str]:
    """Declared scripts that do not exist, sorted."""
    return sorted(
        {s for vs in validators.values() for _, s in vs if not Path(s).is_file()}
    )


def verdict_status(exit_code: Optional[int]) -> str:
    """``pass`` / ``fail``, or ``error`` when the validator itself could not run."""
    if exit_code == 0:
        return "pass"
    if exit_code in (126, 127) or exit_code is None:
        return "error"
    return "fail"


def script_sha256(path: str) -> str:
    """sha256 of a script, ``""`` if it cannot be read."""
    try:
        return hashlib.sha256(Path(path).read_bytes()).hexdigest()
    except OSError:
        return ""


def build_targets(
    produced: Iterable[Tuple[str, str, str]],
    validators: ValidatorMap,
    out_dir: Path,
) -> Dict[str, List[str]]:
    """``{key: [relpath, interp, script, sha256]}``, one entry per (file, validator).

    *produced* is an iterable of ``(stage_id, output_id, file_path)`` for files
    that exist on disk. ``relpath`` is relative to *out_dir* (where snakemake
    runs); ``key`` is ``{relpath}.{i}``, *i* the validator's position on the
    output, and names the result ``.validation/{key}.json``.
    """
    shas = {
        script: script_sha256(script) for vs in validators.values() for _, script in vs
    }
    targets: Dict[str, List[str]] = {}
    for stage_id, output_id, file_path in produced:
        rel = os.path.relpath(str(file_path), str(out_dir))
        for i, (interp, script) in enumerate(validators.get((stage_id, output_id), [])):
            targets[f"{rel}.{i}"] = [rel, interp, script, shas[script]]
    return targets


# Static body of the generated Snakefile. Snakemake placeholders are kept verbatim
# (no str.format here); only the TARGETS literal and __ENV_LINE__ are injected.
_SNAKEFILE_BODY = r'''
rule all:
    input:
        [f".validation/{t}.json" for t in TARGETS]


rule validate:
    input:
        lambda wc: TARGETS[wc.key][0]
    output:
        ".validation/{key}.json"
    log:
        ".validation/{key}.log"
__ENV_LINE__    params:
        interp=lambda wc: TARGETS[wc.key][1],
        script=lambda wc: TARGETS[wc.key][2],
        script_json=lambda wc: json.dumps(TARGETS[wc.key][2]),
        sha=lambda wc: TARGETS[wc.key][3],
        input_json=lambda wc: json.dumps(TARGETS[wc.key][0]),
        log_json=lambda wc: json.dumps(f".validation/{wc.key}.log"),
    shell:
        r"""
        set +e
        if [ -f {params.script:q} ]; then
            {params.interp:q} {params.script:q} {input:q} > {log:q} 2>&1
            code=$?
        else
            echo "validator not found: {params.script}" > {log:q}
            code=127
        fi
        set -e
        printf '{{"output_file": %s, "exit_code": %d, "log": %s, "script": %s, "script_sha256": "%s"}}\n' {params.input_json:q} "$code" {params.log_json:q} {params.script_json:q} {params.sha:q} > {output:q}
        """
'''


def build_snakefile(targets: Dict[str, List[str]], env_directive: str = "") -> str:
    """Render the Snakefile text.

    *env_directive* is a single directive line (e.g. ``conda: "/abs/env.yaml"``)
    or ``""`` to run in the current environment (host backend). It is injected
    into the ``rule validate`` block.
    """
    header = (
        "# Auto-generated validators workflow — do not edit by hand.\n"
        "# Run from the benchmark out/ dir.\n"
        "import json\n"
        f"TARGETS = {json.dumps(targets, indent=4)}\n"
    )
    env_line = f"    {env_directive}\n" if env_directive else ""
    return header + _SNAKEFILE_BODY.replace("__ENV_LINE__", env_line)


def write_validators_snakefile(
    targets: Dict[str, List[str]], env_directive: str, path: Path
) -> None:
    """Write the Snakefile to *path* (creating parents)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(build_snakefile(targets, env_directive))
