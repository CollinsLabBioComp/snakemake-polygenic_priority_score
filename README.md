# Description

This pipeline runs [PoPS](https://github.com/FinucaneLab/pops) analysis for a given set of features. See PoPS GitHub for more information on PoPS analysis.

It performs the following steps:

1. Annotate genes with MAGMA scores
2. Run PoPS with gene features
3. Empirically calculate p-values

## Quickstart

Quickstart for deploying this pipeline locally and on a high performance compute cluster.

### 1. Set up the environment

Create conda environment: `envs/environment.yml` and activate:

```
conda activate snakemake_pps
```

Alternatively, if using singularity or docker, one can pull the image from henryjt/snakemake-polygenic_priority_score:latest.

### 2. Prepare the input files

See `demo/data/features` for example input files. See `docs/README-params.md` for description of parameters.

### 3. Run pipeline

See `demo/run_pops__local_docker.sh` for example of how to run full pipeline

Author: Henry Taylor