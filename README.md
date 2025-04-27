## omnibenchmark

<p align="center">
<a href="https://github.com/omnibenchmark/omnibenchmark/tree/refs/heads/main"><img alt="tests result" src="https://github.com/omnibenchmark/omnibenchmark/workflows/Tests/badge.svg?branch=main"></a>
<a href="https://github.com/omnibenchmark/omnibenchmark"><img alt="Coverage Status" src="./reports/coverage.svg"></a>
<a href="https://github.com/omnibenchmark/omnibenchmark/blob/main/LICENSE"><img alt="License: Apache 2.0" src="https://img.shields.io/badge/License-Apache_2.0-blue.svg"></a>
<a href="https://github.com/psf/black"><img alt="Code style: black" src="https://img.shields.io/badge/code%20style-black-000000.svg"></a>
</p>

[Omnibenchmark](https://omnibenchmark.org), a continuous benchmarking tool.

## Install

With poetry or pip. With micromamba if handling software with conda. See [the tutorial](https://omnibenchmark.org/tutorial/).

### Releases

See [our tags](https://github.com/omnibenchmark/omnibenchmark/tags).

## Developer notes

### omni-schema

Please note the [omni-schema](https://github.com/omnibenchmark/omni-schema) dependency. Benchmark YAML schemas are updated by:

```
- Update src/omni_schema/schema/omni_schema.yaml manually
- `make all`
- Consider `make deploy`
```

And that omni-schema versions are tagged and pinned at the pytoml level ([example](https://github.com/omnibenchmark/omnibenchmark/blob/2ce768bb2cfb693f3e555f751979093964eef63b/pyproject.toml#L38)), so omni-schema changes must precede omnibenchmark changes.

### Documentation

Omnibenchmark docs are served at https://omnibenchmark.org as generated on [Renku's GitLab](https://gitlab.renkulab.io/omnibenchmark/omni_site) with a [review/staging -> production flow](https://gitlab.renkulab.io/omnibenchmark/omni_site/-/blob/master/.gitlab-ci.yml?ref_type=heads).

Documentation includes a CLI reference. This reference is generated/automated [via mkdocs-click](https://gitlab.renkulab.io/omnibenchmark/omni_site/-/blob/master/docs/reference.md?ref_type=heads) and extracts the current CLI commands from [omnibenchmark's `main` head](https://gitlab.renkulab.io/omnibenchmark/omni_site/-/blob/master/requirements.txt?ref_type=heads#L7). Hence, changes to omnibenchmark must be merged to `main` and precede changes to omnibenchmark's documentation.

## Acknowledgements

Omnibenchmark incorporates great FOSS components, including but not limited to: Snakemake, easybuild, apptainer, lmod, LinkML, git. Thank you!

## Preprints

- [Omnibenchmark (alpha) for continuous and open benchmarking in bioinformatics](https://arxiv.org/abs/2409.17038) (2024)
- [Building a continuous benchmarking ecosystem in bioinformatics](https://arxiv.org/abs/2409.15472) (2024)
