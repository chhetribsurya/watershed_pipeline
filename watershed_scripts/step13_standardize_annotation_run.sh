#!/usr/bin/env bash

# GENE LEVEL ANNOTATION
annotation_dir="/scratch16/abattle4/surya/datasets/WatershedAFR/data/annotation/annotation_inputfiles/genelevel_output"
#conda activate r-env

for eachfile in $(ls $annotation_dir/*collapse.tsv);do

    echo -e "\nProcessing Gene level annotation file: $eachfile" 

    # Extract filename without extension
    basefile=$(basename $eachfile .tsv)
    
    # Create the desired output filename
    outfilename="$annotation_dir/${basefile}.stdscaled.tsv"

    Rscript --vanilla ./step13.1_standardize_annotation.R --ANNOT $eachfile --OUTFILE $outfilename

done


# VARIANT LEVEL ANNOTATION
annotation_dir="/scratch16/abattle4/surya/datasets/WatershedAFR/data/annotation/annotation_inputfiles/variantlevel_output"
#conda activate r-env

#for eachfile in $(ls $annotation_dir/*uncollapsed*.tsv);do
for eachfile in $(ls $annotation_dir/gene-GLOBAL-rv.gencode.uncollapsed.rvpair.tsv);do

    echo -e "\nProcessing Variant level annotation file: $eachfile" 

    # Extract filename without extension
    basefile=$(basename $eachfile .tsv)
    
    # Create the desired output filename
    outfilename="$annotation_dir/${basefile}.stdscaled.tsv"

    Rscript --vanilla ./step13.2_standardize_annotation_variantLevel.R --ANNOT $eachfile --OUTFILE $outfilename

done
