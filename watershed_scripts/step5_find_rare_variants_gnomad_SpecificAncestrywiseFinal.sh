#!/bin/bash

# This script performs the end-to-end process of identifying rare variants from gnomAD data,
# specifically for given ancestries. It filters variants, processes them by populations, and
# confirms their rarity against a reference dataset.

# Usage:
# bash step5_find_rare_variants_gnomad_SpecificAncestrywiseFinal.sh -s [scripts_dir] -d [data_dir] -g [kgpex_vcf] -r [regions] -f [gnomad_vcf] -l [pop_list] -p [pop]

# Define command line arguments
while getopts s:d:g:r:f:l:p: flag; do
    case "${flag}" in
        s) scripts_directory=${OPTARG};;
        d) rare_variants_directory=${OPTARG};;
        g) kgpex_raw_vcf=${OPTARG};;
        r) regions_bed=${OPTARG};;
        f) gnomad_raw_vcf=${OPTARG};;
        l) population_list=${OPTARG};;
        p) population=${OPTARG};;
    esac
done

# Check for necessary tools
check_tools() {
    required_tools=("bcftools" "bedtools" "parallel")
    for tool in "${required_tools[@]}"; do
        if ! command -v $tool &> /dev/null; then
            echo "$tool is required but not found. Exiting."
            exit 1
        fi
    done
}

# Function to merge and prepare regions for filtering
merge_regions() {
    regions_merged_name="$(basename $regions_bed .bed)_merged.bed"
    regions_merged=${rare_variants_directory}/${regions_merged_name}
    bedtools merge -i $regions_bed > $regions_merged
    echo $regions_merged
}

# Function to filter variants based on regions and SNP types
filter_variants() {
    local vcf_input=$1
    local vcf_output=$2
    local regions=$3

    if [ ! -f "$vcf_output" ]; then
        bcftools view --regions-file $regions --types snps -Oz -o $vcf_output $vcf_input
        bcftools index --tbi $vcf_output
    fi
}

# Function to subset and recompute allele frequencies for a specific population
subset_and_recompute_af() {
    local vcf_input=$1
    local vcf_output=$2
    local population=$3
    local population_list=$4

    if [ ! -f "$vcf_output" ]; then
        bcftools view --samples-file $population_list $vcf_input | bcftools +fill-tags -Oz -o $vcf_output -- -t AF
        bcftools index --tbi $vcf_output
    fi
}

# Function to filter for rare variants based on allele frequency (AF < 0.01)
filter_rare_variants() {
    local vcf_input=$1
    local vcf_output=$2

    if [ ! -f "$vcf_output" ]; then
        bcftools view --include 'INFO/AF < 0.01' -Oz -o $vcf_output $vcf_input
    fi
}

# Function to confirm rarity of GTEx variants in gnomAD
confirm_rarity_in_gnomad() {
    local kgpex_rare_vcf=$1
    local gnomad_vcf=$2
    local kgpex_rare_bed=$3
    local gnomad_common_vcf=$4
    local gnomad_common_bed=$5
    local kgpex_rare_confirmed_vcf=$6

    if [ ! -f "$kgpex_rare_bed" ]; then
        bcftools query -f '%CHROM\t%POS0\t%END\t%ID\n' $kgpex_rare_vcf > $kgpex_rare_bed
    fi

    if [ ! -f "$gnomad_common_vcf" ]; then
        bash ${scripts_directory}/step5.01_get_gnomad_common_SpecificAncestrywise.sh \
            -g $gnomad_vcf \
            -r $kgpex_rare_bed \
            -d ${rare_variants_directory}/gnomad_by_chr \
            -p $population \
            -o $gnomad_common_vcf
    fi

    if [ ! -f "$gnomad_common_bed" ]; then
        bcftools query -f '%CHROM\t%POS0\t%END\t%ID\n' $gnomad_common_vcf > $gnomad_common_bed
    fi

    if [ ! -f "$kgpex_rare_confirmed_vcf" ]; then
        bedtools intersect -v -a $kgpex_rare_vcf -b $gnomad_common_bed -header | \
        bcftools convert --output $kgpex_rare_confirmed_vcf --output-type z
    fi
}

# Function to get list of rare variants per each gene-individual pair
get_rare_variants_per_gene_individual() {
    local kgpex_rare_confirmed_vcf=$1
    local regions_bed=$2
    local indiv_at_rv=$3
    local rv_sites_raw=$4
    local population=$5
    local rare_variants_directory=$6
    local scripts_directory=$7

    bcftools query -f'[%CHROM\t%POS0\t%END\t%REF\t%ALT\t%INFO/AF\t%SAMPLE\n]' --include 'GT="alt"' $kgpex_rare_confirmed_vcf \
    > $indiv_at_rv

    bedtools intersect -wa -wb -a $indiv_at_rv -b $regions_bed > $rv_sites_raw

    Rscript ${scripts_directory}/step5.02_gene_indiv_rare_variants.R \
    --rv_sites=$rv_sites_raw \
    --popname=$population \
    --outdir=$rare_variants_directory
}

# Main processing sequence
check_tools
regions_merged=$(merge_regions)

kgpex_filtered_vcf=${rare_variants_directory}/kgpex_filtered.vcf.gz
gnomad_filtered_vcf=${rare_variants_directory}/gnomad_filtered.vcf.gz

filter_variants $kgpex_raw_vcf $kgpex_filtered_vcf $regions_merged
filter_variants $gnomad_raw_vcf $gnomad_filtered_vcf $regions_merged

kgpex_pop_vcf=${rare_variants_directory}/kgpex_${population}.vcf.gz
subset_and_recompute_af $kgpex_filtered_vcf $kgpex_pop_vcf $population $population_list

kgpex_pop_rare_vcf=${rare_variants_directory}/kgpex_${population}_rare.vcf.gz
filter_rare_variants $kgpex_pop_vcf $kgpex_pop_rare_vcf

kgpex_rare_bed=${rare_variants_directory}/kgpex_${population}_rare.bed
gnomad_common_vcf=${rare_variants_directory}/gnomad_${population}_common.vcf.gz
gnomad_common_bed=${rare_variants_directory}/gnomad_${population}_common.bed
kgpex_rare_confirmed_vcf=${rare_variants_directory}/kgpex_${population}_rare_confirmed.vcf.gz

confirm_rarity_in_gnomad $kgpex_pop_rare_vcf $gnomad_filtered_vcf $kgpex_rare_bed $gnomad_common_vcf $gnomad_common_bed $kgpex_rare_confirmed_vcf

indiv_at_rv=${rare_variants_directory}/kgpex_${population}_rare_confirmed.indiv.txt
rv_sites_raw=${rare_variants_directory}/gene-${population}-rv.raw.txt

get_rare_variants_per_gene_individual $kgpex_rare_confirmed_vcf $regions_bed $indiv_at_rv $rv_sites_raw $population $rare_variants_directory $scripts_directory

echo "Processing complete. Output available in $rare_variants_directory"

