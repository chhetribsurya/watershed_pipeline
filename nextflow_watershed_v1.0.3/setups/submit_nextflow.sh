#!/usr/bin/env bash

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
  mamba env create -f ../env_ymls/watershed_env.yml
fi

# Activate the environment
conda activate watershed_env


if [ -z "$conda_env_path" ]; then
  echo "Error: Conda environment 'watershed_env' not found."
  echo "'watershed_env' with singularity features required"
  exit 1
fi

# Set environment variables for Nextflow
export NXF_HOME=../cache/nextflowcache
export NXF_WORK=./work
export NXF_TEMP=./tmp
export NXF_CONDA_CACHEDIR=../cache/condacache
export SINGULARITY_CACHEDIR=../cache/singularitycache
#export SINGULARITY_TMPDIR=../cache/singularitycache/tmp
#export NXF_SINGULARITY_LIBRARYDIR=../cache/singularitycache/cachelib

mkdir -p $NXF_HOME
mkdir -p $NXF_WORK
mkdir -p $NXF_TEMP
mkdir -p $NXF_CONDA_CACHEDIR
mkdir -p $SINGULARITY_CACHEDIR

setup_cachedir="$(readlink -f ../cache/setup_cache)"
cachedir_vep="${setup_cachedir}/vepcache"
#conda_cache_dir="${setup_cache_dir}/condacache"
#singularity_cache_dir="${setup_cache_dir}/singularitycache"

mkdir -p $setup_cachedir
mkdir -p $cachedir_vep

#mkdir -p $conda_cache_dir
#mkdir -p $singularity_cache_dir

#nextflow run main.nf --setup_cache_dir "${setup_cachedir}" -profile conda_and_singularity -process.echo -resume
nextflow run main-gnomad.nf --setup_cache_dir "${setup_cachedir}" -profile conda_and_singularity -process.echo -resume

