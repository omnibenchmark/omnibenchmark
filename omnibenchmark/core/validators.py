"""Output validators: discover validator scripts, generate the Snakefile that runs
them against the exact produced files, and read back per-output pass/fail results.

Convention: validators live at ``{validators_dir}/{stage_id}/{output_id}/validate.*``
keyed by the plan's stage id and *output id* (the stable identifiers, e.g.
``PCA/pcas_tsv``), and receive a single produced output file path. They signal
pass/fail by EXIT CODE (0 = pass, non-zero = fail); stdout/stderr are captured to
a per-output log.

Targets are the concrete files the benchmark produced (from the status machinery),
matched to a validator by ``(stage_id, output_id)`` via the output's path template.
We do NOT glob the filesystem: real output trees are deeply nested under hidden
``.{param}`` dirs that globbing won't traverse.

This module is pure (no click / snakemake imports) so it stays easy to test. The
generated Snakefile is run from the benchmark out/ dir; results and logs land under
``out/.validation/`` mirroring each output's path.
"""

import fnmatch
import json
import os
import re
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

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


def path_to_glob(path_template: str) -> str:
    """Turn an output path template into a filename glob.

    ``{dataset}_pcas.tsv`` -> ``*_pcas.tsv``; ``{dataset}.h5ad`` -> ``*.h5ad``.
    """
    return re.sub(r"\{[^}]*\}", "*", path_template)


def build_targets(
    produced: Iterable[Tuple[str, str]],
    validators: ValidatorMap,
    stage_outputs: Dict[str, Dict[str, str]],
    out_dir: Path,
) -> Dict[str, List[str]]:
    """Match produced files to validators -> ``{relpath: [interp, script]}``.

    *produced* is an iterable of ``(stage_id, file_path)`` for files that exist on
    disk. *stage_outputs* maps ``stage_id -> {output_id: path_template}``. A file
    matches a validator when it's under the validator's stage and its basename
    matches the output's filename glob. ``relpath`` is relative to *out_dir* (where
    snakemake runs).
    """
    # Precompute (stage, output_id) -> (filename_glob, spec), only for declared outputs.
    globs: Dict[Tuple[str, str], Tuple[str, Tuple[str, str]]] = {}
    for (stage_id, output_id), spec in validators.items():
        template = stage_outputs.get(stage_id, {}).get(output_id)
        if template:
            globs[(stage_id, output_id)] = (path_to_glob(Path(template).name), spec)

    targets: Dict[str, List[str]] = {}
    for stage_id, file_path in produced:
        base = os.path.basename(file_path)
        for (vstage, _output_id), (glob, spec) in globs.items():
            if vstage == stage_id and fnmatch.fnmatch(base, glob):
                rel = os.path.relpath(str(file_path), str(out_dir))
                targets[rel] = list(spec)
                break
    return targets


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
    shell:
        r"""
        set +e
        {params.interp} {params.script} {input} > {log} 2>&1
        code=$?
        set -e
        printf '{{"output_file": "%s", "exit_code": %d, "log": "%s"}}\n' "{input}" "$code" "{log}" > {output}
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
