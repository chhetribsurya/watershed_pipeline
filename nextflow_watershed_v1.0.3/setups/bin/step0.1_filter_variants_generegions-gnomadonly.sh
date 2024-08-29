#!/bin/bash

# Parse input arguments
while getopts o:r:g: flag; do
    case "${flag}" in
        o) output_dir=${OPTARG};;
        r) gene_regions_file=${OPTARG};;
        g) gnomad_raw_vcf=${OPTARG};;
    esac
done

echo "Setting up environment and verifying tools..."
required_tools=("bcftools" "bedtools" "parallel")
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

# Merge regions for faster filtering
regions_merged="$output_dir/$(basename $gene_regions_file .bed)_merged.bed"
bedtools merge -i $gene_regions_file > $regions_merged
echo "Regions merged into: $regions_merged"

# Filtering Cohort and gnomAD variants
gnomad_vcf="$output_dir/gnomad_filtered.vcf.gz"

# Check and filter gnomAD VCF
if [ ! -f "$gnomad_vcf" ]; then
    bcftools view --regions-file $regions_merged --types snps -Oz -o $gnomad_vcf $gnomad_raw_vcf
    bcftools index --tbi $gnomad_vcf
    echo "gnomAD VCF filtered and indexed at: $gnomad_vcf"
else
    echo "gnomAD VCF already processed."
fi

# Check and filter gnomAD VCF
#if [ ! -f "$gnomad_vcf" ]; then
#    bcftools view --regions-file $regions_merged --types snps -Oz -o $gnomad_vcf $gnomad_raw_vcf
#    bcftools index --tbi $gnomad_vcf
#    echo "gnomAD VCF filtered and indexed at: $gnomad_vcf"
#elif [ "$gnomad_raw_vcf.tbi" -ot "$gnomad_raw_vcf" ]; then  # Check if index is older than VCF
#    echo "Re-indexing gnomAD VCF..."
#    tabix -p vcf $gnomad_raw_vcf
#    echo "gnomAD VCF re-indexed."
#else
#    echo "gnomAD VCF already processed."
#fi


echo "*** Filtering process completed ***"

