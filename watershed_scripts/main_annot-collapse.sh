#!/usr/bin/env bash

# Array of population variables
POPULATIONS=("GLOBAL" "AFR")
POPULATIONS=("GLOBAL" "EUR")
POPULATIONS=("EUR" "AFR")
POPULATIONS=("AFR")
POPULATIONS=("EUR")
POPULATIONS=("GLOBAL" "EUR" "AFR")

# Save log files including rarevariant stats in directory
log_files="./logfiles"
watershed_codedir = "/scratch16/abattle4/surya/datasets/WatershedAFR/WatershedAFR/code/Watershed"

# Check if directory exists
if [ ! -d "$log_files" ]; then
    mkdir -p "$log_files"
fi

# Annotation collapsing for watershed
for POP in "${POPULATIONS[@]}"; do
    echo -e "\n\n***** Annotation collapsing of POP: $POP *****\n"
    #sbatch step12_collapse_annotations_run.sh "$POP"
done


## This part taken care during global outlier calls in newly adapted scripts 
## Collapse/Generate PValue based reference expression Outlier file
##for POP in "${POPULATIONS[@]}"; do
##    echo -e "\n\n***** RV  Expr. outlier reference generation of POP: $POP *****\n"
##    ###sbatch step14_standardize_expOutlierPval_annotation_run.sh ${POP}
##done


# N2pairs annotation generation
for POP in "${POPULATIONS[@]}"; do
    echo -e "\n\n***** N2pairs annotation generation of POP: $POP *****\n"
    sbatch step15_generateN2pairs-final.sh  ${POP}
done


# Merge all collapsed annotations and merge N2pair annotation
for POP in "${POPULATIONS[@]}"; do
    echo -e "\n\n***** Merge all annotations and N2pairs for POP: $POP *****\n"
    sbatch step16_merge_annotations_and_sortN2pair.sh  ${POP}
done

# Watershed prediction
for POP in "${POPULATIONS[@]}"; do
    echo -e "\n\n***** Watershed model prediction for POP: $POP *****\n"
    sbatch ${watershed_codedir}/step17_watershed_model_call.sh  ${POP}
done


