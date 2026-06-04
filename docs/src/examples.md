## Clustering algorithms (*clustbench*)

We have ported the *clustbench* clustering benchmark [(Gagolewski, 2022)](https://clustering-benchmarks.gagolewski.com/index.html) to evaluate **105 datasets** with a known ground truth and **27 methods** using six partition metrics (see table below).

If present, we included noisy points during the clustering process, but ignored them when calculating the performance metrics.

We grouped datasets and methods according to their generator and/or software environment, hence writing modules able to run multiple methods on demand. The [Clustering_conda.yml](https://github.com/omnibenchmark/clustering_example/blob/main/Clustering_conda.yml) manifest makes use of this parametrization to specify the dataset/method/metric to be run from a benchmarking module. 

Beyond *Conda*, we also designed *EasyBuild* and *Apptainer* execution environments to evaluate the impact of the software backend on benchmarking results, both in terms of algorithmic outcomes (e.g., clustering metrics) and computational performance.

### Components

| Stage      | Module        | Components                                                                                                                   | Count |
|------------|---------------|------------------------------------------------------------------------------------------------------------------------------|-------|
| Data       | fcps          | atom, chainlink, engytime, hepta, lsun, target, tetra, twodiamonds, wingnut                                                  | 9     |
|            | graves        | graves1--graves12                                                                                                            | 12    |
|            | other         | aggregation, aniso, blobs, circles, complex9v1--complex9v55                                                                  | 59    |
|            | sipu          | a1, a2, a3, a4, dim032, dim064, dim128, dim256, g2--g6, s1--s4, unbalance, triangle1, triangle2                              | 20    |
|            | uci           | iris, wine, yeast                                                                                                            | 3     |
|            | wut           | spiral, zigzag_outliers                                                                                                      | 2     |
| Clustering | fastcluster   | complete, ward, average, weighted, median, centroid                                                                          | 6     |
|            | sklearn       | birch, kmeans, spectral, gm                                                                                                  | 4     |
|            | agglomerative | average, complete, ward                                                                                                      | 3     |
|            | genieclust    | genie, gic, ica                                                                                                              | 3     |
|            | fcps          | Minimax, MinEnergy, HDBSCAN_2/4/8, Diana, Fanny, Hardcl, Softcl, Clara, PAM                                                  | 11    |
| Metrics    | partition_metrics | normalized_clustering_accuracy, adjusted_fm_score, adjusted_mi_score, adjusted_rand_score, fm_score, mi_score            | 6     |

### Running *clustbench* with Omnibenchmark

To run the benchmark using *Conda* as a software backend:

```bash
git clone git@github.com:omnibenchmark/clustering_example.git
cd clustering_example

ob run Clustering_conda.yml
```

## Community Benchmarks

These are benchmarks organized by the community using Omnibenchmark.

| Git repository | Description | API Version | License | DOI |
|----------------|-------------|-------------|---------|-----|
| [OB_GSEA](https://github.com/cphgeno/OB_GSEA-main) | Comparison of 17 single-sample gene set enrichment analysis (GSEA) methods across diverse datasets | 0.3.2 | MIT | — |
| [scrna-bench](https://github.com/scrna-bench/pipelines-plan) | End-to-end scRNA-seq analysis pipelines, from HDF5 inputs to PCA, clustering, and other outputs | 0.4.0 | MIT | [10.5281/zenodo.19886347](https://zenodo.org/records/19886347), [10.64898/2026.05.01.722166](https://doi.org/10.64898/2026.05.01.722166) |
| [CyTOF pipelines](https://github.com/kaae-2/ob-pipeline-cytof) | Comparison of CyTOF preprocessing and clustering pipelines | 0.3.2 | MIT | - |
| [VarCallBench](https://github.com/Jiayi-Wang-Joey/VarCallBench) | Comparison of variant callers for long read RNA-seq data | 0.5.1 | MIT | [10.5281/zenodo.19857089](https://doi.org/10.5281/zenodo.19857089), [10.64898/2026.04.29.721619](https://doi.org/10.64898/2026.04.29.721619) |
| [OB_LONGREAD_GERMLINE_SV_CALLERS_BENCHMARK](https://github.com/cphgeno/OB_LONGREAD_GERMLINE_SV_CALLERS_BENCHMARK-main) | Benchmarking long-read germline structural variant callers for research and clinic using both linear and pangenome-graph references with a novel dataset | 0.3.2 | MIT | — |
| [OB_SV_BENCHMARK_GERMLINE_SHORTREADS](https://github.com/cphgeno/OB_SV_BENCHMARK_GERMLINE_SHORTREADS-main) | An Omnibenchmark workflow running multiple structural-variant (SV) callers on short-reads data and benchmarking their results against a single truth set with Truvari | 0.3.2 | MIT | — |

Have a benchmark to share? Open a PR (or an [issue](https://github.com/omnibenchmark/omnibenchmark/issues)) to add it to this section.
