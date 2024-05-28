#!/bin/bash

#SBATCH --job-name=pvalN2pairannotation
#SBATCH --output=./logfiles/%x_%j.out
#SBATCH --error=./logfiles/%x_%j.err
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=36G

#POP="EUR"

# Check if an argument is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <POPULATION>"
    exit 1
fi

# Population variable from the argument
POP="$1"

# File paths
datadir1="/scratch16/abattle4/surya/datasets/WatershedAFR/data"
#codedir1="/scratch16/abattle4/surya/datasets/WatershedAFR/WatershedAFR/code/preprocessing/annotation"

rv_file="${datadir1}/rare_variants_gnomad/${POP}/gene-${POP}-rv.txt"
#exp_outlier_file="${datadir1}/data_prep/PEER_${POP}/kgpexV3.${POP}.outlier.controls.v3ciseQTLs_globalOutliersRemoved.txt"
#exp_outlier_file="${datadir1}/rv_expoutlier_refbase/${POP}/MAGE.${POP}.Zscorethresbased.globalOutliersRemoved.txt"
exp_outlier_file="${datadir1}/rv_expoutlier_refbase/${POP}/MAGE.${POP}.Pvalthresbased.globalOutliersRemoved.txt"

#annodir1="${datadir1}/annotation/${POP}/n2pair"
#annodir1="${datadir1}/annotation/${POP}/n2pair_zscore"
annodir1="${datadir1}/annotation/${POP}/n2pair_pval"
output_file="${annodir1}/gene-${POP}-rv.N2pairs.tsv"

# Create n2pair  directory if it doesn't exist
target_dir="${annodir1}"

# Check if directory exists
if [ ! -d "$target_dir" ]; then
    mkdir -p "$target_dir"
fi

# Activate conda env
env_name="/home/schhetr1/anaconda3/envs/py3-env"
eval "$(conda shell.bash hook)"
conda activate "$env_name"

# Run n2pair annotation
echo -e "*** N2pair annotation Processing of Pop: $POP ***"
python ./step15_generateN2pairs-final.py \
    --rare_variant_file $rv_file \
    --exp_outlier_file $exp_outlier_file \
    --outfile $output_file

echo -e "\n *** N2pair annotation completed for Pop: $POP ...\n"
