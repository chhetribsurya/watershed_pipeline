#!/bin/bash

# Parse input arguments
while getopts o:g:r:b:c: flag; do
    case "${flag}" in
        o) output_dir=${OPTARG};;         # Directory where all outputs will be saved
        g) gnomad_common_bed=${OPTARG};;  # Path to the combined gnomAD common bed (converted from VCF)
        r) cohort_rare_vcf=${OPTARG};;    # Path to the cohort rare VCF from previous steps
        b) gene_regions_bed=${OPTARG};;   # BED file of gene regions for mapping
        c) cohort_name=${OPTARG};;        # Cohort name (used for population specificity)
    esac
done

echo "Setting up environment..."
required_tools=("bcftools" "bedtools")
for tool in "${required_tools[@]}"; do
    if ! command -v $tool &> /dev/null; then
        echo "$tool could not be found"
        exit 1
    fi
done
echo "*** Setup and initial checks completed ***"

# Define output files
#gnomad_common_bed="${output_dir}/gnomad_${cohort_name}_common.bed"
cohort_rare_qc_vcf="${output_dir}/cohort_${cohort_name}_rareQC.vcf.gz"
indiv_at_rv="${output_dir}/cohort_${cohort_name}_rareQC.indiv.txt"
rv_sites_raw="${output_dir}/gene-${cohort_name}-rv.raw.txt"

# Check total available cores and determine the number to use
#num_threads=$(nproc)

# Set number of threads
num_threads=$(($(nproc) / 4))
echo "Using $num_threads threads for parallel processing."



# Convert gnomAD common VCF to BED format
#if [ ! -f "$gnomad_common_bed" ]; then
#    echo "Converting gnomAD common VCF to BED..."
#    bcftools query -f '%CHROM\t%POS0\t%END\t%ID\n' $gnomad_common_vcf > $gnomad_common_bed
#    echo "Conversion completed: $gnomad_common_bed"
#else
#    echo "File already exists: $gnomad_common_bed"
#fi


# Select rare variants from cohort based on gnomad common variants match
if [ ! -f "$cohort_rare_qc_vcf" ]; then
    echo "Processing intersection of cohort rare VCF with gnomAD common BED..." 
    #echo -e "\nand bcf tool conversion using the number of threads: $num_threads"
    bedtools intersect -v -a $cohort_rare_vcf -b $gnomad_common_bed -header | \
    bcftools convert --output $cohort_rare_qc_vcf --output-type z
    #bcftools convert --output $cohort_rare_qc_vcf --output-type z --threads $num_threads
    echo "Intersection completed: $cohort_rare_qc_vcf"
else
    echo "File already exists: $cohort_rare_qc_vcf"
fi

# Get list of rare variants per each gene-individual pair
# Position is 0-based like the start position used in bed file format
if [ ! -f "$indiv_at_rv" ]; then
    echo "Getting list of rare variants per each gene-individual pair..."
    echo -e "\n*** List samples that have each rare variant..."
    bcftools query -f'[%CHROM\t%POS0\t%END\t%REF\t%ALT\t%INFO/AF\t%SAMPLE\n]' --include 'GT="alt"' $cohort_rare_qc_vcf > $indiv_at_rv
    echo "List obtained and saved to: $indiv_at_rv"
else
    echo "File already exists: $indiv_at_rv"
fi

# Map rare variants to genes
if [ ! -f "$rv_sites_raw" ]; then
    echo "Mapping rare variants to genes..."
    bedtools intersect -wa -wb -a $indiv_at_rv -b $gene_regions_bed > $rv_sites_raw
    echo "Mapping completed: $rv_sites_raw"
else
    echo "File already exists: $rv_sites_raw"
fi

echo "*** Processing completed ***"

