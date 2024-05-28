#!/usr/bin/env bash

#SBATCH --job-name=global_watershed_run
#SBATCH --output=./logfiles/%x_%j.out
#SBATCH --error=./logfiles/%x_%j.err
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=40G

# Check if an argument is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <POPULATION>"
    exit 1
fi

# Population variable from the argument
POP="$1"
#POP="EUR"

# Activate conda env
env_name="/home/schhetr1/anaconda3/envs/watershed"
eval "$(conda shell.bash hook)"
conda activate "$env_name"

# Run using Watershed approximate inference
# If the number of dimensions (E) is less than equal to 4, exact inference should be used
model="RIVER"  # Can take on "Watershed_exact", "Watershed_approximate", "RIVER"
model="Watershed_approximate"  # Can take on "Watershed_exact", "Watershed_approximate", "RIVER"
model="Watershed_exact"  # Can take on "Watershed_exact", "Watershed_approximate", "RIVER"
number_of_dimensions="3" # Can take on any real number greater than or equal to one
number_of_dimensions="1" # Can take on any real number greater than or equal to one

codedir="/scratch16/abattle4/surya/datasets/WatershedAFR/WatershedAFR/code/Watershed"
datadir="/scratch16/abattle4/surya/datasets/WatershedAFR/data"

#output_prefix="./annotpvalthresh_${POP}_results_new_bonf_stdscaled_annotfilt/${POP}_model_"$model"_number_of_dimensions_"$number_of_dimensions
#input_file="${datadir}/annotation/${POP}/final-${POP}-rv.mergedannotation.plusN2pair.annotfiltered.tsv" # Input file
#input_file="${datadir}/annotation/${POP}/final-${POP}-rv.mergedannotation.plusN2pair.Pvalthresbased.tsv" # Input file
#input_file="${datadir}/annotation/${POP}/final-${POP}-rv.mergedannotation.plusN2pair.signPvalthresbased.stdscaled.tsv" # Input file

# Best 2 methods to test
#input_file="${datadir}/annotation/${POP}/final-${POP}-rv.mergedannotation.plusN2pair.Pvalthresbased.stdscaled.annotfiltered.tsv" # Input file
input_file="${datadir}/annotation/${POP}/final-${POP}-rv.mergedannotation.plusN2pair.Pvalthresbased.stdscaled.tsv" # Input file


# Create n2pair  directory if it doesn't exist
#target_dir="./annotpvalthresh_${POP}_results_new_bonf_stdscaled_annotfilt"
target_dir="./annotpvalthresh_results_new_bonf_stdscaled/${POP}"

# Check if directory exists
if [ ! -d "$target_dir" ]; then
    mkdir -p "$target_dir"
fi

# Set output directory
output_prefix="${target_dir}/${POP}_model_"$model"_number_of_dimensions_"$number_of_dimensions

# RUN WATERSHED EVALUATION
echo -e "\n\n**** WATERSHED EVALUATION STARTS ****"
echo -e "\nWatershed evaluation of: $input_file" 
#Rscript evaluate_watershed.R --input $input_file --number_dimensions $number_of_dimensions --output_prefix $output_prefix --model_name $model --n2_pair_pvalue_fraction 0.009
#Rscript $codedir/evaluate_watershed.R --input $input_file --number_dimensions $number_of_dimensions --output_prefix $output_prefix --model_name $model
#Rscript $codedir/evaluate_watershed.R --input $input_file --number_dimensions $number_of_dimensions --output_prefix $output_prefix --model_name $model --binary_pvalue_threshold 0.0003521127 --n2_pair_pvalue_fraction 0.002115568
#Rscript $codedir/evaluate_watershed.R --input $input_file --number_dimensions $number_of_dimensions --output_prefix $output_prefix --model_name $model
#Rscript $codedir/evaluate_watershed.R --input $input_file --number_dimensions $number_of_dimensions --output_prefix $output_prefix --model_name $model --binary_pvalue_threshold 0.00035
Rscript $codedir/evaluate_watershed.R --input $input_file --number_dimensions $number_of_dimensions --output_prefix $output_prefix --model_name $model --binary_pvalue_threshold 2.535111e-06


#######################
# Run 'predict_watershed.R'
# Note for convenience, the training file is the same as the prediction file. This does not necessarily have to be the case
#######################


# RUN WATERSHED PREDICTION
echo -e "\n\n\n**** WATERSHED PREDICTION STARTS ****"
echo -e "\nWatershed prediction of: $input_file" 
#Rscript $codedir/predict_watershed.R --training_input $input_file --prediction_input $input_file --number_dimensions $number_of_dimensions --output_prefix $output_prefix --model_name $model
#Rscript $codedir/predict_watershed.R --training_input $input_file --prediction_input $input_file --number_dimensions $number_of_dimensions --output_prefix $output_prefix --model_name $model --binary_pvalue_threshold 0.0003521127 
#Rscript $codedir/predict_watershed.R --training_input $input_file --prediction_input $input_file --number_dimensions $number_of_dimensions --output_prefix $output_prefix --model_name $model
Rscript $codedir/predict_watershed.R --training_input $input_file --prediction_input $input_file --number_dimensions $number_of_dimensions --output_prefix $output_prefix --model_name $model --binary_pvalue_threshold 2.535111e-06



