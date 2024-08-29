#!/usr/bin/env bash

# Check if an argument is provided
#if [ $# -ne 1 ]; then
#    echo "Usage: $0 <POPULATION>"
#    exit 1
#fi

# Population variable from the argument
pop="$1"
input_file="$2"
outputdir="$3"
outprefix="$4"
number_of_dim="$5"
model="$6"
watershed_codedir="$7"

# Run using Watershed approximate inference
# If the number of dimensions (E) is less than equal to 4, exact inference should be used
#model="RIVER"  # Can take on "Watershed_exact", "Watershed_approximate", "RIVER"
#model="Watershed_approximate"  # Can take on "Watershed_exact", "Watershed_approximate", "RIVER"
#model="Watershed_exact"  # Can take on "Watershed_exact", "Watershed_approximate", "RIVER"
#number_of_dimensions="3"
#number_of_dimensions="1" 

# Check if directory exists
if [ ! -d "$outputdir" ]; then
    mkdir -p "$outputdir"
fi

# Set output directory
#output_prefix="${output_dir}/${POP}_model_"$model"_number_of_dimensions_"$number_of_dime
output_prefix="${outputdir}/${outprefix}"

# RUN WATERSHED EVALUATION
echo -e "\n\n**** WATERSHED EVALUATION STARTS ****"
echo -e "\nWatershed evaluation of: $input_file" 
Rscript ./evaluate_watershed.R --input $input_file --number_dimensions $number_of_dim --output_prefix $output_prefix --model_name $model --binary_pvalue_threshold 0.01


#-----------------------------------------------------------------------------
# Run 'predict_watershed.R'
# Note for convenience, the training file is the same as the prediction file 
# This does not necessarily have to be the case
#----------------------------------------------------------------------------


# RUN WATERSHED PREDICTION
echo -e "\n\n\n**** WATERSHED PREDICTION STARTS ****"
echo -e "\nWatershed prediction of: $input_file" 
Rscript ./predict_watershed.R --training_input $input_file --prediction_input $input_file --number_dimensions $number_of_dim --output_prefix $output_prefix --model_name $model --binary_pvalue_threshold 0.01

