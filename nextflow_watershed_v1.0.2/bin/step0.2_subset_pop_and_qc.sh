#!/bin/bash

# Parse input arguments
while getopts o:v:c:p: flag; do
    case "${flag}" in
        o) output_dir=${OPTARG};;
        v) cohort_vcf=${OPTARG};;  # This is the filtered cohort VCF from the first script
        c) cohort_name=${OPTARG};;
        p) pop_list=${OPTARG};;
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

# Subset Cohort VCF by population and recompute allele frequencies
#cohort_pop_vcf="$output_dir/${cohort_name}_subset.vcf.gz"
cohort_pop_vcf="$output_dir/$(basename $cohort_vcf .vcf.gz)_subset.vcf.gz"

mkdir -p $output_dir
sample_ids="${output_dir}/${cohort_name}_sample_ids.txt"
awk '{print $1}' $pop_list > $sample_ids

if [ ! -f "$cohort_pop_vcf" ]; then
    #bcftools view --samples-file $pop_list $kgpex | bcftools +fill-tags -Oz -o $kgpex_pop -- -t AF
    bcftools view --samples-file $sample_ids $cohort_vcf | bcftools +fill-tags -Oz -o $cohort_pop_vcf -- -t AF
    bcftools index --tbi $cohort_pop_vcf
    echo "Subsetting Cohort VCF by population completed at: $cohort_pop_vcf"
else
    echo "Subset Cohort VCF by population already exists."
fi

echo "*** Subsetting and QC completed ***"

