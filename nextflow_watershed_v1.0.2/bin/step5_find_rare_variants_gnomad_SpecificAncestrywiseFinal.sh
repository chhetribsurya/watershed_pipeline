#!/bin/bash

# Read in arguments
# while getopts d:g:r:f:p:o: flag
while getopts s:d:g:r:f:l:p:o: flag
do
    case "${flag}" in
        s) scriptdir=${OPTARG};;
        d) rvdir=${OPTARG};;
        g) kgpex_raw=${OPTARG};;
        r) regions=${OPTARG};;
        f) gnomad_raw=${OPTARG};;
        l) pop_list=${OPTARG};;
        p) pop=${OPTARG};;
    esac
done

echo $scriptdir
echo $rvdir
echo $kgpex_raw
echo $regions
echo $gnomad_raw
echo $pop_list;
echo $pop

## Check if bcftools and bedtools are available. If not, exit the script
if ! command -v bcftools &> /dev/null
then
    echo "bcftools could not be found"
    exit
fi

if ! command -v bedtools &> /dev/null
then
    echo "bedtools could not be found"
    exit
fi

if ! command -v parallel &> /dev/null
then
    echo "GNU parallel could not be found"
    exit
fi

# Make output directory
mkdir -p $rvdir
#
## Filter for rare variants within 10kb +/- window around gene body of genes that are protein coding and lincRNA coding
# Save the resulting VCFs to `kgpex.vcf.gz` and `gnomad.padded10kb_PCandlinc_only.vcf.gz` in `${datadir}/rare_variants`

echo "*** Filter for rare variants within 10kb +/- window around gene body of genes that are protein coding and lincRNA coding"

# merge regions that overlap and/or book-end for faster filtering
regions_merged_name="$(basename $regions .bed)_merged.bed"
regions_merged=${rvdir}/${regions_merged_name}
echo bedtools merge -i $regions > $regions_merged


# MAGE variants (keep SNPs only)
kgpex=${rvdir}/kgpex.vcf.gz

if [ -f "$kgpex" ]; then
        echo "**** Filtered GTEx VCF already exists"
else
        echo bcftools view --regions-file $regions_merged --types snps -Oz -o $kgpex $kgpex_raw
        
fi

if [ -f "$kgpex.tbi" ]; then
        echo "**** Filtered GTEx VCF already indexed"
else
        echo "**** Indexing GTEx VCF"
        echo bcftools index --tbi $kgpex
fi


# gnomad variants
gnomad=${rvdir}/gnomad.padded10kb_PCandlinc_only.vcf.gz

if [ -f "$gnomad" ]; then
        echo "**** Filtered gnomAD VCF already exists"
else
        echo bcftools view --regions-file $regions_merged --types snps -Oz -o $gnomad $gnomad_raw 

fi


if [ -f "$gnomad.tbi" ]; then
        echo "**** Filtered gnomAD VCF already indexed"
else
        echo "**** Indexing gnomAD VCF"
        echo bcftools index --tbi $gnomad
fi

echo "*** Done"


## Subset GTEx VCF by population
# also recomputes allele frequencies within that population
echo "*** Subset GTEx VCF by population"

kgpex_pop=${rvdir}/kgpex_${pop}.vcf.gz


#---------------------------------------------------------------------------------------------
# BYpass the following run of recomputing allele frequency as this is global pop based
#cp ${rvdir}/kgpex.vcf.gz $kgpex_pop
#---------------------------------------------------------------------------------------------


if [ -f "$kgpex_pop" ]; then
        echo "**** $kgpex_pop already exists"
else
        #echo cp ${rvdir}/kgpex.vcf.gz $kgpex_pop
        #bcftools +fill-tags recomputes AF after we remove samples
        echo "**** Recomputing $pop AF"
        bcftools view --samples-file $pop_list $kgpex | bcftools +fill-tags -Oz -o $kgpex_pop -- -t AF
fi

if [ -f "$kgpex_pop.tbi" ]; then
        echo "**** Filtered kgpex pop VCF already indexed"
else
        echo "**** Indexing kgpex pop VCF"
        bcftools index --tbi $kgpex_pop
fi

echo "*** Done"


## Filter GTEx VCF for rare variants (MAF < 0.01)
# Save the resulting VCFs as `kgpex_${pop}_rare.vcf.gz` in `${datadir}/rare_variants`
echo "*** Filter KGPEx VCF for rare variants (MAF < 0.01)"

kgpex_pop_rare=${rvdir}/kgpex_${pop}_rare.vcf.gz
if [ -f "$kgpex_pop_rare" ]; then
        echo "**** $kgpex_pop_rare already exists"
else
        echo "performing ancestry operation ..."
        bcftools view --include 'AF<0.01 & AF>0' -Oz -o $kgpex_pop_rare $kgpex_pop
        echo "kgpex rare variant filtering operation completed. check final result file: $kgpex_pop_rare"
fi

echo "*** Done"

## Confirm rarity of GTEx variants in gnomAD
# Quality control to check that rare variants in GTEx are also rare in gnomAD population.
# Save the resulting VCFs to `kgpex_${pop}_rare.QC.vcf.gz` in ${datadir}/rare_variants
echo "*** Confirm rarity of KGPEx variants in gnomAD"

kgpex_rare_bed=${rvdir}/kgpex_${pop}_rare.bed
gnomad_common=${rvdir}/gnomad.${pop}_common.vcf.gz
gnomad_common_bed=${rvdir}/gnomad.${pop}_common.bed
kgpex_pop_rareQC=${rvdir}/kgpex_${pop}_rare.QC.vcf.gz

# Filter gnomAD for variants in GTEx that have AF >= 0.01
if [ -f "$kgpex_rare_bed" ]; then
        echo "**** $kgpex_rare_bed already exists"
else    
        echo "performing vcf to bed conversion"
        bcftools query -f '%CHROM\t%POS0\t%END\t%ID\n' $kgpex_pop_rare > $kgpex_rare_bed
fi

if [ -f "$gnomad_common" ]; then
        echo "**** $gnomad_common already exists"
else

        echo "processing chromosome wise bash script for finding gnomad common variants"
        bash ${scriptdir}/step5.01_get_gnomad_common_SpecificAncestrywise.sh \
        -g $gnomad \
        -r $kgpex_rare_bed \
        -d ${rvdir}/gnomad_by_chr \
        -p $pop \
        -o $gnomad_common

fi

echo "gnomad common variant filtering operation completed. check final result file: $gnomad_common"

if [ -f "$gnomad_common_bed" ]; then
        echo "**** $gnomad_common_bed already exists"
else
        echo "processing conversion of gnomad vcf to gnomad common bed"
        bcftools query -f '%CHROM\t%POS0\t%END\t%ID\n' $gnomad_common > $gnomad_common_bed
fi

# Select rare variants from kgpex based on gnomad common variants match
if [ -f "$kgpex_pop_rareQC" ]; then
        echo "**** $kgpex_pop_rareQC already exists"
else
        echo "processing bedtools intersection of kgpex_pop_rare vcf with gnomad_common_bed"
        bedtools intersect -v -a $kgpex_pop_rare -b $gnomad_common_bed -header | \
        bcftools convert --output $kgpex_pop_rareQC --output-type z
fi

echo "*** Done"
echo "kgpex X gnomad common variant intersection operation completed. check final result file: $kgpex_pop_rareQC"

## Get list of rare variants per each gene-individual pair
# Save the resulting gene-individual-variant files to `gene-${pop}-rv.txt` in `${datadir}/rare_variants`
echo "*** Get list of rare variants per each gene-individual pair"

# Make list of samples with each rare variant. Position is 0-based like the start position in bed file format
echo "**** List samples that have each rare variant"

indiv_at_rv=${rvdir}/kgpex_${pop}_rare.QC.indiv.txt

bcftools query -f'[%CHROM\t%POS0\t%END\t%REF\t%ALT\t%INFO/AF\t%SAMPLE\n]' --include 'GT="alt"' $kgpex_pop_rareQC \
> $indiv_at_rv


# Use bedtools to intersect list of rare variants with protein and lincRNA coding genes padded by 10kb around gene body
# This maps the rare variants to the gene. Two genes can share the same rare variant.
echo "**** Subset rare variants that lie within 10kb around protein coding and lincRNA coding genes"
rv_sites_raw=${rvdir}/gene-${pop}-rv.raw.txt
bedtools intersect -wa -wb -a $indiv_at_rv -b $regions > $rv_sites_raw

echo "Subset rare variants that lie within 10kb around protein coding and lincRNA coding genes operation completed. check final result file: $rv_sites_raw"

#Get list of rare variants per each gene-individual pair
echo "**** Building file with list of rare variants per each gene-individual pair"

Rscript ${scriptdir}/step5.02_gene_indiv_rare_variants.R \
--rv_sites=$rv_sites_raw \
--popname=$pop \
--outdir=$rvdir

echo "*** Done"
echo "Building file with list of rare variants per each gene-individual pair operation completed. check final result file: $rvdir"
echo -e  "\n\n"
echo "******Rare Variant Operation DONE!*******"
echo -e "\n"
date +"%r"
echo -e "\n"
