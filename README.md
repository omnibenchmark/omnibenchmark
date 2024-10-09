## omnibenchmark

<p align="center">
<a href="ttps://github.com/omnibenchmark/omni-py"><img alt="Tests" src="./reports/tests.svg"></a>
<a href="https://github.com/omnibenchmark/omni-py/actions"><img alt="Actions Status" src="https://github.com/omnibenchmark/omni-py/workflows/Tests/badge.svg"></a>
<a href="ttps://github.com/omnibenchmark/omni-py"><img alt="Coverage Status" src="./reports/coverage.svg"></a>
<a href="https://github.com/omnibenchmark/omni-py/blob/main/LICENSE"><img alt="License: Apache 2.0" src="https://img.shields.io/badge/License-Apache_2.0-blue.svg"></a>
<a href="https://github.com/psf/black"><img alt="Code style: black" src="https://img.shields.io/badge/code%20style-black-000000.svg"></a>
</p>


[Omnibenchmark](https://omnibenchmark.org), a continuous benchmarking tool.

## Install

With poetry or pip. With micromamba if handling software with conda. See [the tutorial](https://omnibenchmark.org/tutorial/).

### Releases

See [our tags](https://github.com/omnibenchmark/omnibenchmark/tags).

## Developer notes

Please note the [omni-schema](https://github.com/omnibenchmark/omni-schema) dependency. Benchmark YAML schemas are updated by:

```
- Update src/omni_schema/schema/omni_schema.yaml manually
- `make all`
- Consider `make deploy`
```

## Acknowledgements

Omnibenchmark incorporates great FOSS components, including but not limited to: Snakemake, easybuild, apptainer, lmod, LinkML, git. Thank you!

## Preprints

- [Omnibenchmark (alpha) for continuous and open benchmarking in bioinformatics](https://arxiv.org/abs/2409.17038) (2024)
- [Building a continuous benchmarking ecosystem in bioinformatics](https://arxiv.org/abs/2409.15472) (2024)
  
