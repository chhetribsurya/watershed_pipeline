#!/bin/bash

#SBATCH --job-name=Eur2_ancestry_variants
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=40G
#SBATCH --time=2:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# Activate conda environment
env_name="/home/schhetr1/anaconda3/envs/bcftools"

# Initialize Conda
eval "$(conda shell.bash hook)"
conda activate "$env_name"

# Load R module
module load r/4.2.0

# Set variables
#pop="AFR"
pop="EUR"

rootdir1=/scratch16/abattle4/surya/datasets/WatershedAFR
datadir1=${rootdir1}/data
rawdir1=${rootdir1}/raw_data
codedir1="/scratch16/abattle4/surya/datasets/WatershedAFR/WatershedAFR/code"

# Execute the script
#bash ${codedir1}/preprocessing/rare_variants/step5_find_rare_variants_gnomad_SpecificAncestrywiseFinal-test.sh \
bash ${codedir1}/preprocessing/rare_variants/step5_find_rare_variants_gnomad_SpecificAncestrywiseFinal.sh \
-d ${datadir1}/rare_variants_gnomad \
-g ${rawdir1}/1KG/1KGP_731-samples_all.filtered.phased.vcf.gz \
-r ${datadir1}/data_prep/gencode.v38.GRCh38.genes_padded10kb_PCandlinc_only.bed \
-f /data/abattle4/surya/datasets/for_ashis/gnomAD_v2.1/gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz \
-l /scratch16/abattle4/surya/datasets/WatershedAFR/data/data_prep/${pop}_ids.txt \
-p $pop
