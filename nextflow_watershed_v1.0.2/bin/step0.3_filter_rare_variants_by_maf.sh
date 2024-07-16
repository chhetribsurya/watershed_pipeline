#!/bin/bash

# Parse input arguments
while getopts o:v:c: flag; do
    case "${flag}" in
        o) output_dir=${OPTARG};;   # Output directory for rare variants
        v) cohort_pop_vcf=${OPTARG};;  # This is the subsetted cohort VCF from the subset_and_qc.sh script
        c) cohort_name=${OPTARG};;  # Cohort name used to specify the file
    esac
done

echo "Setting up environment and verifying tools..."
required_tools=("bcftools")
for tool in "${required_tools[@]}"; do
    if ! command -v $tool &> /dev/null; then
        echo "$tool could not be found"
        exit 1
    fi
done
echo "*** Setup and initial checks completed ***"

# Define the output file for rare variants based on MAF
#cohort_rare_vcf="$output_dir/${cohort_name}_rare.vcf.gz"
cohort_rare_vcf="$output_dir/$(basename $cohort_pop_vcf .vcf.gz)_rare.vcf.gz"

if [ -f "$cohort_rare_vcf" ]; then
    echo "**** $cohort_rare_vcf already exists"
else
    echo "Performing filtration for rare variants (MAF < 0.01)..."
    bcftools view --include 'AF<0.01 & AF>0' -Oz -o $cohort_rare_vcf $cohort_pop_vcf
    bcftools index --tbi $cohort_rare_vcf
    echo "$cohort_name rare variant filter file: $cohort_rare_vcf"
fi

echo "*** Rare variant filtering process completed ***"

# Define the output BED file
cohort_rare_bed="${output_dir}/$(basename $cohort_pop_vcf .vcf.gz)_rare.bed"

# Check if the BED file already exists
if [ -f "$cohort_rare_bed" ]; then
    echo "File already exists: $cohort_rare_bed"
else
    echo "Performing VCF to BED conversion..."
    bcftools query -f '%CHROM\t%POS0\t%END\t%ID\n' $cohort_rare_vcf > $cohort_rare_bed
    echo "Conversion completed for VCF $cohort_rare_vcf to $cohort_rare_bed"
fi

echo "*** VCF to BED conversion process completed ***"
