#!/bin/bash
#SBATCH --job-name=piece3
#SBATCH --output=piece3_%A_%a.out
#SBATCH --error=piece3_%A_%a.err
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=80G
#SBATCH --time=01:00:00
#SBATCH --array=0-2

# Array of populations
POPULATIONS=("EUR" "AFR" "GLOBAL")

# Select the population based on the Slurm array task ID
POP=${POPULATIONS[$SLURM_ARRAY_TASK_ID]}

# Define directories
DATA_DIR="/scratch16/abattle4/surya/datasets/WatershedAFR/data/data_prep"
RAREVAR_DIR="/scratch16/abattle4/surya/datasets/WatershedAFR/data/rare_variants_gnomad"
CODEDIR="/scratch16/abattle4/surya/datasets/WatershedAFR/WatershedAFR/code"

# Activate conda env
env_name="/home/schhetr1/anaconda3/envs/r-env"
eval "$(conda shell.bash hook)"
conda activate "$env_name"

# Run the rare variant enrichment analysis R script
Rscript ${CODEDIR}/rare_var_enrichment.R \
        ${DATA_DIR} \
        ${POP} \
        ${DATA_DIR}/PEER_${POP}/kgpexV3.${POP}.outlier.controls.v3ciseQTLs_globalOutliersRemoved.txt \
        ${RAREVAR_DIR}/gene-${POP}-rv.txt

