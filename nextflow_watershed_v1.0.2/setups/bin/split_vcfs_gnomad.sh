#!/bin/bash

# Parse input arguments
while getopts o:v: flag; do
    case "${flag}" in
        o) output_dir=${OPTARG};;        # Directory where the chromosome-specific VCFs will be saved
        v) gnomad_filtered_vcf=${OPTARG};; # Path to the filtered gnomAD VCF file
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

# Create output directory
mkdir -p $output_dir
echo "Output directory created at: $output_dir"

# Extract basename without extension for naming output files
base_name=$(basename $gnomad_filtered_vcf .vcf.gz)

# Split the VCF into separate files by chromosome
echo "Splitting VCF by chromosome..."
for chr in {1..22} X Y; do
    output_vcf="${output_dir}/${base_name}_chr${chr}.vcf.gz"
    if [ ! -f "$output_vcf" ]; then
        echo "Processing chromosome: chr$chr..."
        bcftools view -r chr$chr -Oz -o $output_vcf $gnomad_filtered_vcf
        bcftools index -t $output_vcf
        echo "Chromosome chr$chr VCF generated at: $output_vcf"
    else
        echo "Chromosome chr$chr VCF already exists: $output_vcf"
    fi
done

echo "*** VCF splitting process completed ***"

