"""Output validators: discover validator scripts, generate the Snakefile that runs
them against the exact produced files, and read back per-output pass/fail results.

A validator is declared on a stage output (``validator: {script}``) or dropped at
``{validators_dir}/{stage_id}/{output_id}/validate.*``; a declared one wins. It
receives a single produced output file path and signals pass/fail by EXIT CODE
(0 = pass, 126/127 = could not run, anything else = fail); stdout/stderr are
captured to a per-output log.

Targets are the concrete files the benchmark produced (from the status machinery),
keyed by the output id they were produced for, so no filename matching is needed.
We do NOT glob the filesystem: real output trees are deeply nested under hidden
``.{param}`` dirs that globbing won't traverse.

Pure helpers; ``omnibenchmark.backend.validators`` runs them. The generated
Snakefile is run from the benchmark out/ dir; results and logs land under
``out/.validation/`` mirroring each output's path.
"""

import hashlib
import json
import os
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

# Interpreter by validator extension. Anything else falls back to bash.
INTERP = {".py": "python3", ".R": "Rscript", ".sh": "bash"}

VALIDATION_DIR = ".validation"

# Key: (stage_id, output_id). Value: (interpreter, absolute script path).
ValidatorMap = Dict[Tuple[str, str], Tuple[str, str]]


def discover_validators(validators_dir: Path) -> ValidatorMap:
    """Map ``(stage_id, output_id) -> (interpreter, abs_script)``.

    Walks ``{validators_dir}/{stage_id}/{output_id}/validate.*``. Raises if a
    directory holds more than one ``validate.*`` (multiple-per-output is not
    supported yet).
    """
    validators: ValidatorMap = {}
    if not validators_dir.is_dir():
        return validators
    for stage_dir in sorted(p for p in validators_dir.iterdir() if p.is_dir()):
        for out_dir in sorted(p for p in stage_dir.iterdir() if p.is_dir()):
            scripts = sorted(p for p in out_dir.glob("validate.*") if p.is_file())
            if not scripts:
                continue
            if len(scripts) > 1:
                raise ValueError(
                    f"multiple validators in {out_dir}: {[s.name for s in scripts]}"
                )
            script = scripts[0]
            validators[(stage_dir.name, out_dir.name)] = (
                INTERP.get(script.suffix, "bash"),
                str(script.resolve()),
            )
    return validators


def collect_validators(
    model, benchmark_dir: Path, validators_dir: Path
) -> Tuple[ValidatorMap, Set[Tuple[str, str]]]:
    """Discovered validators overlaid with the plan's declared ones.

    Returns ``(validators, unmatched)``; *unmatched* are discovered
    ``(stage, output)`` keys the plan does not declare, usually a typo.
    """
    found = discover_validators(validators_dir)
    plan_keys = {(st.id, o.id) for st in model.stages for o in st.outputs}
    for st in model.stages:
        for o in st.outputs:
            if o.validator:
                script = (benchmark_dir / o.validator.script).resolve()
                found[(st.id, o.id)] = (INTERP.get(script.suffix, "bash"), str(script))
    unmatched = set(found) - plan_keys
    return {k: v for k, v in found.items() if k in plan_keys}, unmatched


def verdict_status(exit_code: Optional[int]) -> str:
    """``pass`` / ``fail``, or ``error`` when the validator itself could not run."""
    if exit_code == 0:
        return "pass"
    if exit_code in (126, 127) or exit_code is None:
        return "error"
    return "fail"


def _sha256(path: str) -> str:
    try:
        return hashlib.sha256(Path(path).read_bytes()).hexdigest()
    except OSError:
        return ""


def build_targets(
    produced: Iterable[Tuple[str, str, str]],
    validators: ValidatorMap,
    out_dir: Path,
) -> Dict[str, List[str]]:
    """Map produced files to validators -> ``{relpath: [interp, script, sha256]}``.

    *produced* is an iterable of ``(stage_id, output_id, file_path)`` for files
    that exist on disk. ``relpath`` is relative to *out_dir* (where snakemake runs).
    """
    shas = {script: _sha256(script) for _, script in validators.values()}
    return {
        os.path.relpath(str(file_path), str(out_dir)): [
            *validators[(stage_id, output_id)],
            shas[validators[(stage_id, output_id)][1]],
        ]
        for stage_id, output_id, file_path in produced
        if (stage_id, output_id) in validators
    }


# Static body of the generated Snakefile. Snakemake placeholders are kept verbatim
# (no str.format here); only the TARGETS literal and __ENV_LINE__ are injected.
_SNAKEFILE_BODY = r'''
wildcard_constraints:
    relpath = r"(?!\.validation/).+"


rule all:
    input:
        [f".validation/{t}.json" for t in TARGETS]


rule validate:
    input:
        "{relpath}"
    output:
        ".validation/{relpath}.json"
    log:
        ".validation/{relpath}.log"
__ENV_LINE__    params:
        interp=lambda wc: TARGETS[wc.relpath][0],
        script=lambda wc: TARGETS[wc.relpath][1],
        script_json=lambda wc: json.dumps(TARGETS[wc.relpath][1]),
        sha=lambda wc: TARGETS[wc.relpath][2],
        input_json=lambda wc: json.dumps(wc.relpath),
        log_json=lambda wc: json.dumps(f".validation/{wc.relpath}.log"),
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
