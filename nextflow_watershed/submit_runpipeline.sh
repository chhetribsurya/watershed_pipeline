#!/bin/bash

# run_pipeline.sh
# This script sets environment variables and then runs the Nextflow pipeline.

# Module load environments
ml java/19

# Retrieve the Conda environment path for 'watershed_pyenv'
conda_env_path=$(conda info --envs | grep watershed_pyenv | awk '{print $2}')

if [ -z "$conda_env_path" ]; then
  echo "Error: Conda environment 'watershed_pyenv' not found."
  exit 1
fi

# Set environment variables for Nextflow
export NXF_HOME=./cache/nextflowcache
export NXF_WORK=./work
export NXF_TEMP=./tmp
export NXF_CONDA_CACHEDIR=./cache/condacache
export NXF_IGNORE_RESUME_HISTORY=True

# Create directories if they don't exist
mkdir -p "$NXF_HOME"
mkdir -p "$NXF_WORK"
mkdir -p "$NXF_TEMP"
mkdir -p "$NXF_CONDA_CACHEDIR"

chmod -R 755 "$SINGULARITY_CACHEDIR"

# Absolute path to the bin directory
BIN_DIR=$(pwd)/bin

# Path for conda YAML file
CONDA_PYENV_YML=$(pwd)/env_ymls/watershed_pyenv.yml
CONDA_RENV_YML=$(pwd)/env_ymls/watershed_renv.yml
CONDA_ENV_YML=$(pwd)/env_ymls/watershed_env.yml

# Create cache dir
CACHE_DIR=$(pwd)/cache

# Run the Nextflow pipeline
nextflow run ./main.nf -c nextflow.config \
  --analysis STRATIFIED \
  --bin_dir "$BIN_DIR" \
  --cache_dir "$CACHE_DIR" \
  --watershed_pyenv "$CONDA_PYENV_YML" \
  --watershed_renv "$CONDA_RENV_YML" \
  --watershed_env "$CONDA_ENV_YML" \
  --genotype_pcs 5 \
  -profile conda \
  -process.echo \
  -resume

# Run the Nextflow pipeline with correct Conda environment and a specific run name for resuming
nextflow run ./main.nf -c nextflow.config \
  --analysis GLOBAL \
  --cohort_name GLOBAL \
  --bin_dir "$BIN_DIR" \
  --cache_dir "$CACHE_DIR" \
  --watershed_pyenv "$CONDA_PYENV_YML" \
  --watershed_renv "$CONDA_RENV_YML" \
  --watershed_env "$CONDA_ENV_YML" \
  --genotype_pcs 5 \
  -profile conda \
  -process.echo \
  -resume

#---------------------------------------
