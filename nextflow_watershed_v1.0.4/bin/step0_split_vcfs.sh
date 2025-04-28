#!/bin/bash

# Usage: ./split_vcfs.sh -v /path/to/input.vcf.gz -d /path/to/output_directory

while getopts v:d: flag; do
    case "${flag}" in
        v) vcf_path=${OPTARG};;   # Path to the input VCF file
        d) output_dir=${OPTARG};; # Directory where chromosome VCFs will be saved
    esac
done

echo "Setting up environment..."
required_tools=("bcftools")
for tool in "${required_tools[@]}"; do
    if ! command -v $tool &> /dev/null; then
        echo "$tool could not be found"
        exit 1
    fi
done
echo "*** Setup and initial checks completed ***"

# Create output directory if it doesn't exist
mkdir -p $output_dir
echo "Output directory created at: $output_dir"

# Extract the base name of the VCF file
#base_name=$(basename "$vcf_path" .vcf.gz)

# Edited version
base_name=$(basename "$vcf_path")
base_name=${base_name%.vcf.gz}
base_name=${base_name%.vcf.bgz}

# Extract the list of chromosomes from the VCF file
chromosomes=$(bcftools index --stats $vcf_path | cut -f1 | uniq)

echo "Starting split of VCF by chromosome..."

# Loop over each chromosome and extract it into a separate VCF file
for chr in $chromosomes; do
    output_vcf="${output_dir}/${base_name}_${chr}.vcf.gz"
    if [ -f "$output_vcf" ]; then
        echo "File already exists: $output_vcf"
    else
        echo "Processing chromosome: $chr"
        bcftools view -r $chr -O z -o $output_vcf $vcf_path
        bcftools index --tbi $output_vcf
        echo "Chromosome $chr extracted to $output_vcf"
    fi
done

echo "*** Chromosome-wise VCF split completed ***"

