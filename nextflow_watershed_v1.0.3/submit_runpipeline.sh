#!/bin/bash

# This script sets all environment variables and runs the Nextflow pipeline

# Module load environments
ml java/19
ml mamba
ml anaconda/2020.07 

# Load Conda commands
eval "$(conda shell.bash hook)"

conda_env_path=$(conda env list | grep watershed_env | awk '{print $1}'| head -n 1)

# Check if the environment name/path has been captured
if [ -n "$conda_env_path" ]; then
  echo "Environment 'watershed_env' already exists, activating it..."
else
  echo "Creating new environment 'watershed_env'..."
  mamba env create -f ./env_ymls/watershed_env.yml
fi

# Activate the env
conda activate watershed_env


if [ -z "$conda_env_path" ]; then
  echo "Error: Conda environment 'watershed_env' not found."
  echo "'watershed_env' with singularity features required"
  exit 1
fi

# Set environment variables for pipeline
export NXF_HOME=./cache/nextflowcache
export NXF_WORK=./work
export NXF_TEMP=./tmp
export SINGULARITY_CACHEDIR=./cache/singularitycache
#export NXF_SINGULARITY_LIBRARYDIR=./cache/singularitycache/cachelib
#export NXF_CONDA_CACHEDIR=./cache/condacache
export NXF_IGNORE_RESUME_HISTORY=True

# Create directories if they don't exist
mkdir -p "$NXF_HOME"
mkdir -p "$NXF_WORK"
mkdir -p "$NXF_TEMP"
mkdir -p "$SINGULARITY_CACHEDIR"
#mkdir -p "$NXF_SINGULARITY_LIBRARYDIR"
#mkdir -p "$NXF_CONDA_CACHEDIR"
chmod -R 755 "$SINGULARITY_CACHEDIR"

# Absolute path to the bin directory
BIN_DIR=$(pwd)/bin

# Path for conda YAML file
CONDA_PYENV_YML=$(pwd)/env_ymls/watershed_pyenv.yml
CONDA_RENV_YML=$(pwd)/env_ymls/watershed_renv.yml
CONDA_CADDENV_YML=$(pwd)/env_ymls/watershed_caddenv.yml
CONDA_ENV_YML=$(pwd)/env_ymls/watershed_env.yml

# Create cache dir
CACHE_DIR=$(pwd)/cache

SETUP_CACHEDIR="$(readlink -f ./cache/setup_cache)"

SETUP_CACHEDIR="$CACHE_DIR/setup_cache"
VEP_CACHEDIR="${SETUP_CACHEDIR}/vepcache"

## Running the Nextflow Pipeline
nextflow run ./main.nf -c nextflow.config \
  --cohort_name gtexEUR \
  --ancestry EUR \
  --bin_dir "$BIN_DIR" \
  --cache_dir "$CACHE_DIR" \
  --watershed_pyenv "$CONDA_PYENV_YML" \
  --watershed_renv "$CONDA_RENV_YML" \
  --watershed_caddenv "$CONDA_CADDENV_YML" \
  --watershed_env "$CONDA_ENV_YML" \
  --genotype_pcs 5 \
  --skip_cache \
  --cachedir_vep "$VEP_CACHEDIR" \
  --tpm_infile "/scratch16/abattle4/surya/datasets/test_vep/test_watershed_scripts/gtex_input/tpm_Adipose-Subcutaneous.tab" \
  --readcount_infile "/scratch16/abattle4/surya/datasets/test_vep/test_watershed_scripts/gtex_input/read_count_Adipose-Subcutaneous.tab" \
  --covariate_infile "/scratch16/abattle4/surya/datasets/test_vep/test_watershed_scripts/gtex_input/covariates_Adipose-Subcutaneous.tab" \
  --subjids_file "/scratch16/abattle4/surya/datasets/test_vep/test_watershed_scripts/gtex_input/subjids_Adipose-Subcutaneous_EUR.txt" \
  --rv_file "/scratch16/abattle4/surya/datasets/test_vep/test_watershed_scripts/gtex_input/Adipose-Subcutaneous_subsample.vcf.gz" \
  --tissue "Adipose" \
  -profile conda \
  -resume

  #--rv_file "/scratch16/abattle4/surya/datasets/project_watershed/nextflow_pipe/gtex_runs/gtex_input/Adipose-Subcutaneous.vcf.bgz" \

