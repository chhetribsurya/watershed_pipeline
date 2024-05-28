#!/bin/bash

#SBATCH --job-name=GencodeAnnot
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=20G

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

gencode_file="/scratch16/abattle4/surya/datasets/WatershedAFR/raw_data/1KG/gencode.v38.GRCh38.genes.gtf"
rv_file="${datadir1}/rare_variants_gnomad/${POP}/gene-${POP}-rv.txt"
annodir1="${datadir1}/annotation/${POP}/gencode"
output_file="${annodir1}/gene-${POP}-rv.gencode.txt"

# Create gencode annotation directory if it doesn't exist
if [ ! -d "${annodir1}" ]; then
    mkdir -p "${annodir1}"
fi

# Activate conda env
env_name="/home/schhetr1/anaconda3/envs/py3-env"
eval "$(conda shell.bash hook)"
conda activate "$env_name"

# Run Gencode annotation to compute dist to TSS and TES from rarevars
echo -e "*** Gencode Annot Processing of Pop: $POP ***"
echo -e "\n\nScriptrun: python ./step10_gencode_tss_tes_distance.py $gencode_file $rv_file $output_file \n\n"
python ./step10_gencode_tss_tes_distance.py "$gencode_file" "$rv_file" "$output_file"

