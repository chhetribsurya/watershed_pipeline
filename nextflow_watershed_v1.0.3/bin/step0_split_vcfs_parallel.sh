#!/bin/bash

# Usage: ./split_vcfs_parallel.sh -v /path/to/input.vcf.gz -d /path/to/output_directory

#while getopts v:d: flag; do
while getopts v:d:c: flag; do
    case "${flag}" in
        v) vcf_path=${OPTARG};;   # Path to the input VCF file
        d) output_dir=${OPTARG};; # Directory where chromosome VCFs will be saved
        c) num_of_cores=${OPTARG};; # Number of Cores
    esac
done

#num_of_cores=${num_of_cores:-8} # Default to 8 cores if not set

# Automatically determine the number of available CPU cores
#num_of_cores=$(nproc)
#echo "Detected $num_of_cores cores. Using all available cores..."

# Check total available cores and determine the number to use
total_cores=$(nproc)
#num_of_cores=${user_cores:-$total_cores}
num_of_cores=$(( num_of_cores > total_cores ? total_cores : num_of_cores ))

echo "Using $num_of_cores cores for parallel processing..."


echo "Setting up environment..."
required_tools=("bcftools" "parallel")
for tool in "${required_tools[@]}"; do
    if ! command -v $tool &> /dev/null; then
        echo "$tool could not be found"
        exit 1
    fi
done

echo "*** Setup and initial checks completed ***"

# Check if the VCF file is indexed, index if not
if [ ! -f "${vcf_path}.csi" ] && [ ! -f "${vcf_path}.tbi" ]; then
    echo "Indexing VCF file..."
    bcftools index $vcf_path
    echo "Indexing completed."
else
    echo "VCF file is already indexed."
fi

# Create output directory if it doesn't exist
mkdir -p $output_dir
echo "Output directory created at: $output_dir"

# Extract the base name of the VCF file
base_name=$(basename "$vcf_path" .vcf.gz)

# List all chromosomes from the VCF file and store in an array
mapfile -t chromosomes < <(bcftools index --stats $vcf_path | cut -f1 | uniq)

# Display the number of cores to be used for parallel processing
echo "Using $num_of_cores cores for parallel processing..."

# Use parallel to process each chromosome
echo "Starting split of VCF by chromosome using parallel..."

#parallel --will-cite -j "$num_of_cores" \
#"bcftools view -r {} -Oz -o ${output_dir}/${base_name}_{}.vcf.gz ${vcf_path} && \
#bcftools index --tbi ${output_dir}/${base_name}_{}.vcf.gz" ::: ${chromosomes[@]}

#parallel --will-cite -j $num_of_cores \
#    "echo 'Processing chromosome: {}' && \
#     bcftools view -r {} -Oz -o ${output_dir}/${base_name}_{}.vcf.gz ${vcf_path} && \
#     bcftools index --tbi ${output_dir}/${base_name}_{}.vcf.gz && \
#     echo 'Completed chromosome: {}'" ::: ${chromosomes[@]}

parallel --will-cite -j $num_of_cores \
    "output_vcf=${output_dir}/${base_name}_{}.vcf.gz; \
     if [[ -f \$output_vcf ]]; then \
        echo 'File already exists: \$output_vcf'; \
     else \
        echo 'Processing chromosome: {}'; \
        bcftools view -r {} -Oz -o \$output_vcf ${vcf_path} && \
        bcftools index --tbi \$output_vcf && \
        echo 'Completed chromosome: {}'; \
     fi" ::: ${chromosomes[@]}

echo -e "\n*** Chromosome-wise VCF split completed ***\n"


