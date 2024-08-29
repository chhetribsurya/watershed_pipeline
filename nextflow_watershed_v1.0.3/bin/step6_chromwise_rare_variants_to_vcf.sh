#!/bin/bash

# CONVERT RARE VARIANT TO VCF chromosomewise for CADD run

# Check if an argument is provided
if [ $# -ne 4 ]; then
    echo "Usage: $0 <POPULATION> <RVFILE> <OUTPUTDIR> <SCRIPTDIR>"
    exit 1
fi

# Population variable from the argument
POP="$1"
rarevar_file="$2"
chromwise_dir="$3"
script_dir="$4"

#rarevar_file="${rarevar_dir}/gene-${POP}-rv.txt"

# Check if the $POP directory exists. If not, create it.
#rarevar_pop_dir=${rarevar_dir}/${POP}
#if [ ! -d "${rarevar_pop_dir}" ]; then
#    mkdir "${rarevar_pop_dir}"
#fi

# Move rarevar_pop file and reset rarevar file directory
# mv ${rarevar_file} ${rarevar_pop_dir}/gene-${POP}-rv.txt
# rarevar_file="${rarevar_pop_dir}/gene-${POP}-rv.txt"

# Set chromosome-wise directory based on $POP name 
# rarevar_pop_dir=$(dirname "${rarevar_file}")
# chromwise_dir="${rarevar_pop_dir}/chromwise"

# Check if the chromwise directory exists. If not, create it.
if [ ! -d "${chromwise_dir}" ]; then
    mkdir "${chromwise_dir}"
fi

# Check for a header and save it to a variable
header=$(head -n 1 "${rarevar_file}")

# Loop through the unique chromosomes in the rarevar_file
for chrom in $(awk 'NR > 1 {print $3}' "${rarevar_file}" | sort | uniq); do
    # Define output file format
    outfile="${chromwise_dir}/gene-${POP}-rv.${chrom}.txt"
    
    # Print the header into the respective output file
    echo "$header" > "${outfile}"
    
    # Filter rows for the specific chromosome and append to the output file
    awk -v chrom="$chrom" '$3 == chrom {print}' "${rarevar_file}" >> "${outfile}"
    
    # Now, run the R script for each chromosome-specific file
    echo -e "*** Processing file : ${outfile}\n"

    # Run Rscript with output file from above
    Rscript "${script_dir}/step6.1_rare_variants_to_vcf.R" --RV "${outfile}"

done
