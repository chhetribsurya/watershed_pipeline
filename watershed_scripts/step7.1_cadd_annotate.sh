#!/bin/bash

# Load the conda environment
cadd="/home/schhetr1/anaconda3/envs/cadd"
eval "$(conda shell.bash hook)"
conda activate $cadd

# Define the location of CADD scripts
cadd_loc="/scratch16/abattle4/surya/tools/CADD-scripts"

# We'll expect chrom, file, corenum, and pop as arguments
corenum=$1
chrom=$2
pop=$3
file=$4
outdir=$5

# Run CADD command
echo "Executed Command line: bash ${cadd_loc}/CADD.sh -a -g GRCh38 -c ${corenum} -o ${outdir}/gene-${pop}-rv.CADD.${chrom}.tsv.gz ${file}"

bash ${cadd_loc}/CADD.sh -a -g GRCh38 -c ${corenum} -o ${outdir}/gene-${pop}-rv.CADD.${chrom}.tsv.gz ${file}

# Deactivate the conda environment
conda deactivate
