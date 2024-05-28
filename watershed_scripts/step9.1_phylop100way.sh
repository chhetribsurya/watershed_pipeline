#!/bin/bash

# Obtain conservation scores for rare variants from PhyloP 100way scores in bedgraph format
## read in arguments
while getopts b:r:p: flag
do
    case "${flag}" in
        b) phylop=${OPTARG};; #/data/annotation/hg38.phyloP100way.bedGraph
        r) rarevariantsbed=${OPTARG};; #/data/rare_variants_gnomad/gene-GLOBAL-rv.bed
        p) pop=${OPTARG};;
    esac
done

# Load the conda environment
watershed="/home/schhetr1/anaconda3/envs/watershed"
eval "$(conda shell.bash hook)"
conda activate $watershed

# defining directories
codedir1="/scratch16/abattle4/surya/datasets/WatershedAFR/WatershedAFR/code/preprocessing/annotation"
annodir=$(dirname $phylop)
outdir=$(dirname $rarevariantsbed)

combinedbed=${annodir}/hg38.phyloP100way.sorted.bed
bgzippedcombined=${combinedbed}.gz
if [ -f "$bgzippedcombined" ]; then

    echo "**** hg38.phyloP100way bgzipped file  already exists"

else
    # split by chromosome
    for i in {1..22};do
          chrom=chr$i
          bed=${annodir}/hg38.phyloP100way.${chrom}.bed
          bedsorted=${annodir}/hg38.phyloP100way.${chrom}.sorted.bed
         echo "Extracting $chrom"
         grep -w $chrom $phylop > $bed

         # sort on each split bed file
         echo "Sorting $chrom"
         sort --parallel=8 -T $annodir -k1,1 -k2,2n $bed > $bedsorted
    done

    # combine by chromosomes
    sortedlist=${annodir}/sorted.bed.files.list
    combinedbed=${annodir}/hg38.phyloP100way.sorted.bed
    echo "Combining sorted bed files by chromosome"
    ls ${annodir}/hg38.phyloP100way.*.sorted.bed | sort -V > $sortedlist
    cat $(grep -v '^#' $sortedlist) > $combinedbed

    # compress with bigzip
    bgzippedcombined=${combinedbed}.gz
    echo "bigzipping"
    bgzip -c $combinedbed > $bgzippedcombined

fi

if [ -f "$bgzippedcombined.tbi" ]; then

    echo "**** hg38.phyloP100way indexed file already exists"

else
    # index with tabix
    echo "indexing"
    tabix -p bed $bgzippedcombined
fi

# query phylop scores for the rare variants
echo "querying phylop 100 way scores for population: $pop"
python ${codedir1}/step9.2_phylop100way_query.py --phylop $bgzippedcombined --rarevariant $rarevariantsbed --pop $pop

