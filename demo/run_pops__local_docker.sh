#!/bin/sh

# Check for docker
docker --version

# Get repo to pass to Docker
SNK_REPO=$(pwd)/../.

# Run command inside Docker image -- snakemake only uses singularity natively
mkdir -p results # make results dir to mount 
docker run --name snakemake-polygenic_priority_score \
     --mount type=bind,source="$SNK_REPO",target=/home/container_user/snakemake-polygenic_priority_score \
     --mount type=bind,source="$SNK_REPO/assets",target=/home/container_user/assets \
     --mount type=bind,source="${PWD}",target=/home/container_user/analysis \
     --rm -t \
     henryjt/snakemake-polygenic_priority_score \
     /bin/bash -c "cd /home/container_user/analysis &&
          source activate snakemake_pps &&
          snakemake \
               --snakefile /home/container_user/snakemake-polygenic_priority_score/Snakefile \
               --configfile config_analysis.yml \
               --printshellcmds \
               --cores 1 \
               --use-conda \
               $1
          "
