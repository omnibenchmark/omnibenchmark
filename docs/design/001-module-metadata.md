# 001: Module Metadata

[![Status: Implemented](https://img.shields.io/badge/Status-Implemented-brightgreen.svg)](https://github.com/omnibenchmark/docs/design)
[![Version: 0.1](https://img.shields.io/badge/Version-0.1-blue.svg)](https://github.com/omnibenchmark/docs/design)

**Authors**: ben
**Date**: 2025-06-18
**Status**: Implemented
**Version**: 0.1.0
**Supersedes**: N/A
**Reviewed-by**: dincicau, csoneson, imallona
**Related Issues**: #145

## Changes

| Version | Date | Description | Author |
|---------|------|-------------|--------|
| 0.1.0 | 2025-06-18 | Initial draft | ben |
| 0.1.0 | 2025-12-16 | Marked as implemented | ben |


## 1. Problem statement

Previous [spec for module metadata](https://github.com/omnibenchmark/internal_docs/blob/master/architecture/method_contributor_design.md) suggested to expose attribution fields (for citation) and execution details (entrypoints, environment specs) in a single YAML file.

We also rely on a `config.cfg` file to be [present at the root of a valid omnibenchmark module](https://github.com/imallona/clustbench_fcps/blob/main/config.cfg). As of today, this file only contains information about the module entrypoint:

```
[DEFAULT]
SCRIPT=run_fcps.R
```

### Non-Goals

Tangentially related, it has also been noted in the past that the `1:1` relation between module and entrypoint leads to a perhaps unnecessary proliferation of repos. I'm noting in here that allowing more than one entry point could allow to reuse a repo for different benchmark stages, since it's just a namespace that can be selected in the benchmark yaml. But objections were raised to "breaking the atomicity paradigm", so we'll not consider the multiple entrypoints in scope for this design document.


### Design Goals

- Leverage existing tooling (GitHub, Zenodo) for citation.
- Lower the barrier to entry for new contributors.
- Keep execution metadata concise, clearly separated.
- Support extensions (new entrypoints, env types) without bloating either file.

## 2. File Roles

Current spec expose both citation (attribution) metadata and execution details (entrypoints, environment specs) in a single YAML. We propose moving to a hybrid spec:

1. `omnibenchmark.yaml` for module execution semantics (entrypoints, named stages, environment specs: conda, easyconfig).

```YAML
entrypoints:
  # just a default entrypoint is permitted for the time being
  default: shuffler.run
environments:
  conda: this-environment.yaml
  easybuild: https://…/AMGX-2.4.0.eb
```

2. Adopt a `CITATION.cff` alongside it for all citation/attribution fields (title, authors, abstract, references, license).

```YAML
cff-version: 1.2.0
message: "If you use this module in your publication, please cite as below."
title: "Shuffler"
abstract: |
  Shuffles (per row) a rectangular gz-compressed, tab-separated text file.
  Known to break on embedded quotes. Reruns with the same seed are identical.
authors:
  - family-names: Turtle
    given-names: John
date-released: 2025-06-13
version: 1.0.0
license: GPL-3.0-or-later
references:
  - title: “Original code”
    authors:
      - family-names: Bird
        given-names: John
    url: http://stolensnippets.com/shuffler
  - title: “Patch and data”
    authors:
      - family-names: Fish
        given-names: John
    doi: 10.1038/s41586-024-07042-7
```


It remains to be defined what omnibenchmark execution should do with the provided environments. One conservative option is to use the benchmark YAML environments, but make use of the contributor provided environment (if any) as references of known, good working environments (in terms of version compatibility, it's understood that module contributors only provide assurances about the tested list of explicit and implicit dependencies.)

### Semantic Fields

| File | Purpose | Key Fields |
|------|---------|------------|
| CITATION.cff | Attribution metadata | title, abstract, authors, date-released, version, license, references |
| omnibenchmark.yaml | Module execution & environment specification | entrypoints (default & named), environments (conda, easyconfig, etc.) |

## 3. Analysis & Alternatives

After evaluating multiple approaches to module metadata, we've identified three main alternatives with distinct trade-offs. The hybrid approach offers the best balance between standardization and flexibility.

### All-in-one YAML
- **Pros**:
  - Single source.
- **Cons**:
  - Duplicates CFF logic.
  - Ignores ecosystem tools.
  - Harder to extend citation features.

### Pure CITATION.cff
- **Pros**:
  - Standard citation.
- **Cons**:
  - Lacks entrypoint/env semantics.
  - No tool support for runtime keys.
  - Risks misusing CFF for execution details.

### Hybrid (proposed)
- **Pros**:
  - Clear separation of concerns.
  - Each file plays to its strengths.
  - Extensible via custom keys in omnibenchmark.yaml.
  - Minimal learning curve.
- **Cons**:
  - Two files, slight overhead.
  - Requires documenting the split.

## 4. Implementation details

In omnibenchmark, one or more commands should:

- Recursively descend all contributed modules
- Validate `CITATION.cff` and `omnibenchmark.yaml`:

    - Must include required execution keys (e.g., default entrypoint)
    - Must contain valid free software license (can be a warning, validate SPDX as we're doing now)
- Generate a joint citation file for the benchmark (concatenated .cff, for Zenodo export or ob cite)

## 5. Extensibility & Future Work

- `omnibenchmark.yaml` can adopt new environments types as they are needed (e.g. Docker, Singularity).
- `omnibenchmark.yaml` can define new entrypoints as they are needed.


## 6. Recommendations and next steps

- Adopt the hybrid approach.
- Make a bare minimum of fields mandatory (entrypoint, URL, author, LICENSE)
- Update docs to guide module contributors, including best practices for generating citation metadata (e.g. using online tools).
- Explicitely address differences, if any found relevant after further analysis, between metadata about code, algorithms, datasets and publications. Separating attribution for distinct named entities might be relevant to automate, e.g. download or discovery of datasets or code (just to mention a plausible use case.)
- Provide templates for both files in module scaffolding generation.
- Consider adding a tool to guide module contributers to generate CITATION.cff via a text-wizard.
- Watch the ecosystem in case new CFF spec are released that we want to incorporate.

## 7. References

1. [Citation File Format (CFF) Specification](https://github.com/citation-file-format/citation-file-format)
2. [Previous Module Metadata Specification](https://github.com/omnibenchmark/internal_docs/blob/master/architecture/method_contributor_design.md)
3. [Zenodo Metadata Documentation](https://developers.zenodo.org/#representation)
4. [Example Omnibenchmark Module: clustbench_fcps](https://github.com/imallona/clustbench_fcps)
