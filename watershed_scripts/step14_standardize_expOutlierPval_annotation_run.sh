#!/usr/bin/env bash

# load R
ml r/4.2.0

# collapse expression outlierfile as a preprocessing first before Pvalue conversion
Rscript step14.0_collapse_expoutlierfile.R

# Column containing Zscore in the annotation file
Zcol="Zscore"

# GENE LEVEL ANNOTATION
annotation_dir="/scratch16/abattle4/surya/datasets/WatershedAFR/data/annotation/annotation_inputFinal/genelevel_output"
#conda activate r-env

for eachfile in $(ls $annotation_dir/*globalOutliersRemoved*collapse*.tsv);do
#for eachfile in $(ls $annotation_dir/kgpexV3.GLOBAL.outlier.controls.v3ciseQTLs_globalOutliersRemoved.collapse.fullgenes.tsv);do

    echo -e "\nProcessing Gene level annotation file: $eachfile" 

    # Extract filename without extension
    basefile=$(basename $eachfile .tsv)
    
    # Create the desired output filename
    outfilename="$annotation_dir/${basefile}.PVAL.tsv"

    Rscript ./step14.1_convert_z_to_p.R --ZFILE $eachfile --OUTFILE $outfilename --ZCOL $Zcol

done


# VARIANT LEVEL ANNOTATION
annotation_dir="/scratch16/abattle4/surya/datasets/WatershedAFR/data/annotation/annotation_inputFinal/variantlevel_output"
#conda activate r-env

for eachfile in $(ls $annotation_dir/*globalOutliersRemoved*uncollapsed*.tsv);do

    echo -e "\nProcessing Variant level annotation file: $eachfile" 

    # Extract filename without extension
    basefile=$(basename $eachfile .tsv)
    
    # Create the desired output filename
    outfilename="$annotation_dir/${basefile}.PVAL.tsv"

    Rscript ./step14.2_convert_z_to_p_variantLevel.R --ZFILE $eachfile --OUTFILE $outfilename --ZCOL $Zcol

done
