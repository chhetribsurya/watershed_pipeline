#!/bin/bash

# This script finds common variants in gnomAD for multiple ancestries.

# Usage:
# bash step5.01_get_gnomad_common_SpecificAncestrywise.sh -g [gnomad_vcf] -r [kgpex_rare_bed] -d [output_dir] -p [pop] -o [output_vcf]

# Define command line arguments
while getopts g:r:d:p:o: flag; do
    case "${flag}" in
        g) gnomad_vcf=${OPTARG};;
        r) kgpex_rare_bed=${OPTARG};;
        d) output_directory=${OPTARG};;
        p) population=${OPTARG};;
        o) output_vcf=${OPTARG};;
    esac
done

# Function to split gnomAD VCF by chromosome
split_gnomad_by_chromosome() {
    local gnomad_vcf=$1
    local output_directory=$2

    mkdir -p $output_directory
    for i in {1..22}; do
        bcftools view -r chr$i --output-file ${output_directory}/gnomad_chr${i}.vcf.gz --output-type z $gnomad_vcf
        bcftools index --tbi ${output_directory}/gnomad_chr${i}.vcf.gz
    done
}

# Function to find common variants in gnomAD by chromosome
find_common_variants_by_chr() {
    local gnomad_chr_vcf=$1
    local kgpex_rare_bed=$2
    local output_directory=$3
    local population=$4

    chr=$(basename ${gnomad_chr_vcf##*_} .vcf.gz)
    population_info=""

    case $population in
        AFR) population_info="INFO/AF_afr >= 0.01";;
        EUR) population_info="INFO/AF_nfe >= 0.01";;
        AMR) population_info="INFO/AF_amr >= 0.01";;
        EAS) population_info="INFO/AF_eas >= 0.01";;
        SAS) population_info="INFO/AF_oth >= 0.01";;
        GLOBAL) population_info="INFO/AF_afr >= 0.01 || INFO/AF_amr >= 0.01 || INFO/AF_eas >= 0.01 || INFO/AF_nfe >= 0.01 || INFO/AF_oth >= 0.01";;
        *) echo "Population is neither AFR, EUR, AMR, EAS, SAS, nor GLOBAL. Exiting."; exit 1;;
    esac

    bcftools view --regions-file $kgpex_rare_bed --include "$population_info" \
        -Oz -o ${output_directory}/gnomad_${population}_common_$chr.vcf.gz $gnomad_chr_vcf
}

# Function to combine chromosome VCFs into a single file
combine_chr_vcfs() {
    local output_directory=$1
    local population=$2
    local output_vcf=$3

    by_chr_list=${output_directory}/by_chr_list.txt
    ls ${output_directory}/gnomad_${population}_common_chr*.vcf.gz | sort -V > $by_chr_list
    bcftools concat --file-list $by_chr_list --output $output_vcf --output-type z
}

# Main processing sequence
split_gnomad_by_chromosome $gnomad_vcf $output_directory

export -f find_common_variants_by_chr
export kgpex_rare_bed
export output_directory
export population

parallel --jobs 12 find_common_variants_by_chr {} $kgpex_rare_bed $output_directory $population ::: ${output_directory}/gnomad_chr*.vcf.gz
combine_chr_vcfs $output_directory $population $output_vcf

echo "Done"
