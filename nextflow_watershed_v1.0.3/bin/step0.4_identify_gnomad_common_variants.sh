#!/bin/bash

# Parse input arguments
while getopts o:g:r:c: flag; do
    case "${flag}" in
        o) output_dir=${OPTARG};;       # Directory to store output VCF and BED files
        g) gnomad_chr_vcf=${OPTARG};;   # Path to a specific gnomAD VCF (already split by chromosome)
        r) cohort_rare_bed=${OPTARG};;  # Path to cohort rare BED from previous step
        c) ancestry=${OPTARG};;         # Ancestry, used for filtering specific populations from gnomAD data
    esac
done

echo "Setting up environment and processing a specific chromosome..."
required_tools=("bcftools")
for tool in "${required_tools[@]}"; do
    if ! command -v $tool &> /dev/null; then
        echo "$tool could not be found"
        exit 1
    fi
done
echo "*** Setup and initial checks completed ***"

# Extract basename without extension
basename=$(basename $gnomad_chr_vcf .vcf.gz)

# Create directory for output files
mkdir -p $output_dir
echo "Output directory prepared at: $output_dir"

# Define output file names using the basename of the input file
gnomad_common_vcf="${output_dir}/${basename}_common.vcf.gz"
gnomad_common_bed="${output_dir}/${basename}_common.bed"
#num_threads=$(nproc)

# Set number of threads
#num_threads=$(($(nproc) / 7))
#echo "Using $num_threads threads for parallel processing."


#gnomad_common_sort_bed="${output_dir}/${basename}_common.srt.bed"

# Filter for common variants based on cohort
echo "Filtering for common variants in gnomAD for file: $basename"
case "$ancestry" in
    "AFR") filter_expr='INFO/AF_afr >= 0.01' ;;
    "EUR") filter_expr='INFO/AF_nfe >= 0.01' ;;
    "AMR") filter_expr='INFO/AF_amr >= 0.01' ;;
    "EAS") filter_expr='INFO/AF_eas >= 0.01' ;;
    "SAS") filter_expr='INFO/AF_oth >= 0.01' ;;
    "GLOBAL") filter_expr='INFO/AF_afr >= 0.01 || INFO/AF_amr >= 0.01 || INFO/AF_eas >= 0.01 || INFO/AF_nfe >= 0.01 || INFO/AF_oth >= 0.01' ;;
    *) echo "Cohort name is neither AFR, EUR, AMR, EAS, SAS, GLOBAL. Exiting"; exit 1 ;;
esac

#bcftools view --regions-file $cohort_rare_bed --include "$filter_expr" -Oz -o $gnomad_common_vcf $gnomad_chr_vcf --threads $num_threads
bcftools view --regions-file $cohort_rare_bed --include "$filter_expr" -Oz -o $gnomad_common_vcf $gnomad_chr_vcf
bcftools index --tbi $gnomad_common_vcf
#bcftools index --tbi $gnomad_common_vcf --threads $num_threads

echo "Filtered VCF prepared at: $gnomad_common_vcf"

# Convert filtered VCF to BED format
if [ ! -f "$gnomad_common_bed" ]; then
    echo "Converting filtered VCF to BED format..."
    bcftools query -f '%CHROM\t%POS0\t%END\t%ID\n' $gnomad_common_vcf > $gnomad_common_bed
    echo "Conversion to BED completed: $gnomad_common_bed"

    # Sort the BED file
    #echo "Sorting the BED file..."
    #bedtools sort -i $gnomad_common_bed > $gnomad_common_sort_bed
    #echo "Sorting completed: $gnomad_common_sort_bed"

else
    echo "BED file already exists: $gnomad_common_bed"
fi

echo "Processing complete for file: $basename"

