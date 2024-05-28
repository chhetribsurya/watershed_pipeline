#!/usr/bin/env bash

#SBATCH --job-name=merge_all_annotations_withN2pair_pvalthresh
#SBATCH --output=./logfiles/%x_%j.out
#SBATCH --error=./logfiles/%x_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=15G
#SBATCH --partition=shared

# Set variables for the job
#POP="GLOBAL"  # Population identifier used in file naming

# Check if an argument is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <POPULATION>"
    exit 1
fi

# Population variable from the argument
POP="$1"

DATADIR="/scratch16/abattle4/surya/datasets/WatershedAFR/data"
ANNOTATION_DIR="${DATADIR}/annotation"  # Directory containing annotation files

#REFBASE_FILE="${DATADIR}/rv_expoutlier_baseref/${POP}/gene-${POP}-rv.expoutlier.Pval.reference.tsv"  # Path to reference base file
#REFBASE_FILE="${DATADIR}/rv_expoutlier_refbase/${POP}/gene-${POP}-rv.expoutlier.Zscorethresbased.reference.tsv"  # Path to reference base file
REFBASE_FILE="${DATADIR}/rv_expoutlier_refbase/${POP}/gene-${POP}-rv.expoutlier.Pvalthresbased.reference.tsv"  # Path to reference base file

#N2PAIR_FILE="${DATADIR}/annotation/${POP}/n2pair/gene-${POP}-rv.N2pairs.tsv"  # Path to N2 pairs file
#N2PAIR_FILE="${DATADIR}/annotation/${POP}/n2pair_zscore/gene-${POP}-rv.N2pairs.tsv"  # Path to N2 pairs file
N2PAIR_FILE="${DATADIR}/annotation/${POP}/n2pair_pval/gene-${POP}-rv.N2pairs.tsv"  # Path to N2 pairs file

# Activate conda env
env_name="/home/schhetr1/anaconda3/envs/py3-env"
eval "$(conda shell.bash hook)"
conda activate "$env_name"

# Merge baseref annotations with n2pair 
mkdir -p $ANNOTATION_DIR
echo -e "*** Merge annotations with N2pair of Pop: $POP ***"
python ./step16_merge_annotations_and_sortN2pair.py \
    --pop "$POP" \
    --annotation_dir "$ANNOTATION_DIR" \
    --refbase_file "$REFBASE_FILE" \
    --n2pair_file "$N2PAIR_FILE"

echo -e "\n *** N2pair annotation merging completed for Pop: $POP ...\n"

# SignPVAL based: Merge baseref annotations with n2pair 
echo -e "*** SignPVAL based: Merge annotations with N2pair of Pop: $POP ***"
python ./step16_merge_annotations_and_sortN2pair_SignPVAL.py \
    --pop "$POP" \
    --annotation_dir "$ANNOTATION_DIR" \
    --refbase_file "$REFBASE_FILE" \
    --n2pair_file "$N2PAIR_FILE"

echo -e "\n *** SignPVAL based: N2pair annotation merging completed for Pop: $POP ...\n"


