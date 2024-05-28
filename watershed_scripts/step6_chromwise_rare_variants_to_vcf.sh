#!/bin/bash

# Directories and other setup parameters
codedir1="/scratch16/abattle4/surya/datasets/WatershedAFR/WatershedAFR/code"
rarevar_file="/scratch16/abattle4/surya/datasets/WatershedAFR/data/rare_variants_gnomad/gene-GLOBAL-rv.txt"
pop="GLOBAL"

# Extract directory from the rarevar_file path
parent_dir=$(dirname "${rarevar_file}")
chromwise_dir="${parent_dir}/chromwise"

# Check if the chromwise directory exists. If not, create it.
if [ ! -d "${chromwise_dir}" ]; then
    mkdir "${chromwise_dir}"
fi

# Check for a header and save it to a variable
header=$(head -n 1 "${rarevar_file}")

# Load R module for Rscript
ml r/4.2.0

# Loop through the unique chromosomes in the rarevar_file
for chrom in $(awk 'NR > 1 {print $3}' "${rarevar_file}" | sort | uniq); do
    
    # Define output file format
    outfile="${chromwise_dir}/gene-${pop}-rv.${chrom}.txt"
    
    # Print the header into the respective output file
    echo "$header" > "${outfile}"
    
    # Filter rows for the specific chromosome and append to the output file
    awk -v chrom="$chrom" '$3 == chrom {print}' "${rarevar_file}" >> "${outfile}"
    
    # Now, run the R script for each chromosome-specific file
    echo -e "*** Processing file : ${outfile}\n"

    # Run Rscript with output file from above
    Rscript "${codedir1}/preprocessing/rare_variants/step6.1_rare_variants_to_vcf.R" --RV "${outfile}"

done
