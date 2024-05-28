#!/bin/bash

# help function
function usage {
    echo "Usage: $0"
    echo "Performs the following steps:"
    echo "1) Activates a conda environment."
    echo "2) bgzip's a VCF file while retaining the original."
    echo "3) Indexes the gzipped VCF with tabix."
    echo "4) Runs a Python script using the gzipped VCF as input."
    exit 1
}

# Check if user asks for help
if [[ "$1" == "--help" || "$1" == "-h" ]]; then
    usage
fi

# Step 1: Activate conda environment
eval "$(conda shell.bash hook)"
conda activate watershed
if [[ $? -ne 0 ]]; then
    echo "Failed to activate watershed conda environment. Exiting."
    exit 1
fi

# Step 2: bgzip the VCF file
input_file="/scratch16/abattle4/surya/datasets/WatershedAFR/data/annotation/vep/gene-GLOBAL-rv.vep.loftee.vcf"

echo -e "\nBlock zipping file: $input_file"
bgzip -c $input_file > "${input_file}.gz"
echo -e "\nBgzip Done ..."

if [[ $? -ne 0 ]]; then
    echo "Failed to bgzip the VCF file. Exiting."
    exit 1
fi

# Step 3: Index the gzipped VCF using tabix
echo -e "\nTabixing file: $input_file"
tabix -p vcf "${input_file}.gz"
echo -e "\nTabix Done ..."

if [[ $? -ne 0 ]]; then
    echo "Failed to index the VCF file using tabix. Exiting."
    exit 1
fi

# Step 4: Run the Python script
vep_loftee_file="${input_file}.gz"
echo -e "\nPython parsing of vep vcf file ..."
python step11.1_parse_vep_loftee_Fullchrom.py --anno $vep_loftee_file

if [[ $? -ne 0 ]]; then
    echo "Failed to run the Python script. Exiting."
    exit 1
fi

echo "All steps completed successfully!"
