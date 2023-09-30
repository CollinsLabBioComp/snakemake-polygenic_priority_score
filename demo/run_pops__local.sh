#!/bin/sh

## Activate QTL conda environment
# source activate gwas_snakemake

## Set any necessary env variables
export BCFTOOLS_PLUGINS="$(conda info --base)/envs/gwas_snakemake/libexec/bcftools"
export SNK_DIR="/Users/taylorhj/repo/snakemake-polygenic_priority_score"

## ADD ANY SPECIFIC ENV NEEDS

snakemake \
     --snakefile ${SNK_DIR}/Snakefile \
     --configfile config_analysis.yml \
     --printshellcmds \
     --cores 1 \
     $1
