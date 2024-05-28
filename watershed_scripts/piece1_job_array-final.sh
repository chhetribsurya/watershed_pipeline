#!/bin/bash

#SBATCH --job-name=piece1_processing
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=20G
#SBATCH --time=01:00:00
#SBATCH --output=./logfiles/%x_%A_%a.out
#SBATCH --error=./logfiles/%x_%A_%a.err
#SBATCH --array=0-2

# Array of populations
POPULATIONS=("EUR" "AFR" "GLOBAL")

# Select the population based on the Slurm array task ID
POP=${POPULATIONS[$SLURM_ARRAY_TASK_ID]}

# Define directories
rootdir1=/scratch16/abattle4/surya/datasets/WatershedAFR
datadir1=${rootdir1}/data
codedir1="/scratch16/abattle4/surya/datasets/WatershedAFR/WatershedAFR/code"

# Activate conda env
env_name="/home/schhetr1/anaconda3/envs/py3-env"
eval "$(conda shell.bash hook)"
conda activate "$env_name"

# PRE-Processing data for outlier call
python step1_process_data-final.py \
        $DATA_DIR \
        $DATA_DIR/${POP}_ids.txt \
        $POP

# Activate conda env:
env_name="/home/schhetr1/anaconda3/envs/r-env"
eval "$(conda shell.bash hook)"
conda activate "$env_name"

# Log2 TPM processing and filtering
Rscript ${codedir1}/preprocessing/data_prep/step1_preprocess_expr.R \
  --COV ${datadir1}/data_prep/PEER_${POP}/kgpex_v3_eQTL_covariates.txt \
  --PEER ${datadir1}/data_prep/PEER_${POP}

# Compute PEER factors
bash ${codedir1}/preprocessing/data_prep/step2.0_calculate_PEER_factors.sh \
-p ${datadir1}/data_prep/PEER_${POP} \
-t ${codedir1}/PEER/peer_cache \
-d ${codedir1}/PEER/peer_1.3.sif

# Compute residuals (does not need PEER package)
bash ${codedir1}/preprocessing/data_prep/step2.1_calculate_PEER_residuals.sh \
  -p ${datadir1}/data_prep/PEER_${POP} \
  -c ${datadir1}/data_prep/kgpex_v3_eQTL_covariates.txt

# Define directories
DATA_DIR="/scratch16/abattle4/surya/datasets/WatershedAFR/data/data_prep"
RAREVAR_DIR="/scratch16/abattle4/surya/datasets/WatershedAFR/data/rare_variants_gnomad"
OUTPUT_DIR="/scratch16/abattle4/surya/datasets/WatershedAFR/data/rv_expoutlier_refbase"

# Expression outlier calling plus global outlier removal
echo -e "\nCMD: Rscript step2_outlier_calling_script-final.R \
    ${DATA_DIR}/PEER_${POP}/LCL_Blood.peer.v3ciseQTLs.ztrans.txt \
    ${RAREVAR_DIR}/${POP}/gene-${POP}-rv.txt \
    ${OUTPUT_DIR}/${POP} \
    $POP \n"

Rscript step2_outlier_calling_script-final.R \
    ${DATA_DIR}/PEER_${POP}/LCL_Blood.peer.v3ciseQTLs.ztrans.txt \
    ${RAREVAR_DIR}/${POP}/gene-${POP}-rv.txt \
    ${OUTPUT_DIR}/${POP} \
    $POP

echo -e "\n\nTASK COMPLETED ...\n\n"


