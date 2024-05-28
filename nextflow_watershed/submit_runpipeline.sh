#!/bin/bash

# run_pipeline.sh
# This script sets environment variables to prevent Nextflow from writing to your home directory
# and then runs the Nextflow pipeline.

# Module load environments
ml java/19
#ml singularity/3.8.7

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
#export SINGULARITY_CACHEDIR=./cache/singularitycache
#export SINGULARITY_TMPDIR=./cache/singularitycache/tmp
#export NXF_SINGULARITY_LIBRARYDIR=./cache/singularitycache/cachelib
export NXF_CONDA_CACHEDIR=./cache/condacache
export NXF_IGNORE_RESUME_HISTORY=True

#Environment Variables Desriptions:
#NXF_HOME: Directory for Nextflow's home files, preventing it from using the default location in your home directory.
#NXF_WORK: Directory for storing intermediate pipeline results.
#NXF_TEMP: Directory for temporary files.
#SINGULARITY_CACHEDIR: Directory for Singularity cache files.

# Create directories if they don't exist
mkdir -p "$NXF_HOME"
mkdir -p "$NXF_WORK"
mkdir -p "$NXF_TEMP"
#mkdir -p "$SINGULARITY_CACHEDIR"
#mkdir -p "$SINGULARITY_TMPDIR"
#mkdir -p "$NXF_SINGULARITY_LIBRARYDIR"
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
#nextflow run ./main.nf -c nextflow.config --bin_dir $BIN_DIR --watershed_pyenv "$conda_env_path" -profile conda 
#nextflow run ./main.nf -c nextflow.config --bin_dir $BIN_DIR --watershed_pyenv "$conda_env_yml" -profile conda -process.echo -bg 
#nextflow run ./main.nf -c nextflow.config --cohort GLOBAL --bin_dir "$BIN_DIR" --cache_dir "$CACHE_DIR" --watershed_pyenv "$CONDA_ENV_YML" -profile conda -process.echo
#nextflow run ./main.nf -c nextflow.config --cohort STRATIFIED --bin_dir "$BIN_DIR" --cache_dir "$CACHE_DIR" --watershed_pyenv "$CONDA_ENV_YML" -profile conda -process.echo
#nextflow run ./main.nf -c nextflow.config --analysis STRATIFIED --bin_dir "$BIN_DIR" --cache_dir "$CACHE_DIR" --watershed_pyenv "$CONDA_ENV_YML" -profile conda

# Final Nextflow Runs
#nextflow run ./main.nf -c nextflow.config --analysis STRATIFIED --bin_dir "$BIN_DIR" --cache_dir "$CACHE_DIR" --watershed_pyenv "$CONDA_PYENV_YML" --watershed_renv "$CONDA_RENV_YML" -profile conda -resume
#nextflow run ./main.nf -c nextflow.config --analysis GLOBAL --cohort_name GLOBALpop --bin_dir "$BIN_DIR" --cache_dir "$CACHE_DIR" --watershed_pyenv "$CONDA_ENV_YML" -profile conda

#---------------------------------------

# Run the Nextflow pipeline with a specific run name for resuming
#-w /scratch16/abattle4/surya/datasets/for_scripts/project_watershed/nextflow_pipe/test-scripts/work_folder \
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
  #-resume -bg > log.txt


#-profile conda,singularity \
#  -w /scratch16/abattle4/surya/datasets/for_scripts/project_watershed/nextflow_pipe/test-scripts/work_folder \

# Run the Nextflow pipeline with correct Conda environment and a specific run name for resuming
#nextflow run ./main.nf -c nextflow.config \
#  --analysis GLOBAL \
#  --cohort_name GLOBAL \
#  --bin_dir "$BIN_DIR" \
#  --cache_dir "$CACHE_DIR" \
#  --watershed_pyenv "$CONDA_PYENV_YML" \
#  --watershed_renv "$CONDA_RENV_YML" \
#  --watershed_env "$CONDA_ENV_YML" \
#  --genotype_pcs 5 \
#  -profile conda \
#  -process.echo \
#  -resume

#---------------------------------------



#  -name stratified_run
#  -name global_run


# Run Nextflow pipeline with the absolute path
#nextflow run ./main.nf \
#    --out_dir ./ \
#    --subjids_file /scratch16/abattle4/surya/datasets/WatershedAFR/data/data_prep/EUR_ids.txt \
#    --covariate_infile /scratch16/abattle4/surya/datasets/WatershedAFR/data/data_prep/eQTL.covariates.tab.gz \
#    --pseudocount_infile /scratch16/abattle4/surya/datasets/WatershedAFR/data/data_prep/kgpex.sample.unfiltered.pseudocounts.tab \
#    --tpm_infile /scratch16/abattle4/surya/datasets/WatershedAFR/data/data_prep/kgpex.sample.unfiltered.tpm.tab \
#    --variable EUR \
#    --tissue LCL \
#    -with-script $BIN_DIR/step1_process_data.py

