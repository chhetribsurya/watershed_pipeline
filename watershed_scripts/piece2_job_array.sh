#!/bin/bash

#SBATCH --job-name=piece2_ancestry_rare_variants
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=40G
#SBATCH --time=22:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --array=0-2

# Array of populations
POPULATIONS=("AFR" "EUR" "GLOBAL")

# Select the population based on the Slurm array task ID
pop=${POPULATIONS[$SLURM_ARRAY_TASK_ID]}

# Set directories
rootdir1="/scratch16/abattle4/surya/datasets/WatershedAFR"
datadir1="${rootdir1}/data"
rawdir1="${rootdir1}/raw_data"
codedir1="${rootdir1}/WatershedAFR/code"

# Activate conda environment
env_name="/home/schhetr1/anaconda3/envs/bcftools"
eval "$(conda shell.bash hook)"
conda activate "$env_name"

# Load R module
module load r/4.2.0

# Execute the script for finding rare variants
bash "${codedir1}/preprocessing/rare_variants/step5_find_rare_variants_gnomad_SpecificAncestrywiseFinal.sh" \
    -d "${datadir1}/rare_variants_gnomad" \
    -g "${rawdir1}/1KG/1KGP_731-samples_all.filtered.phased.vcf.gz" \
    -r "${datadir1}/data_prep/gencode.v38.GRCh38.genes_padded10kb_PCandlinc_only.bed" \
    -f "/data/abattle4/surya/datasets/for_ashis/gnomAD_v2.1/gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz" \
    -l "${datadir1}/data_prep/${pop}_ids.txt" \
    -p "$pop"

