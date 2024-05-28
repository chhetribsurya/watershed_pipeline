#!/usr/bin/env bash

##SBATCH --job-name=N2pairannotation
##SBATCH --output=./logfiles/%x_%j.out
##SBATCH --error=./logfiles/%x_%j.err
##SBATCH --time=2:00:00
##SBATCH --ntasks=1
##SBATCH --cpus-per-task=10
##SBATCH --mem=36G

#POP="EUR"

# Check if an argument is provided
#if [ $# -ne 1 ]; then
#    echo "Usage: $0 <POPULATION>"
#    exit 1
#fi

# Population variable from the argument
#POP="$1"
POP="EUR"

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
input_file="${datadir}/annotation/${POP}/final-${POP}-rv.mergedannotation.plusN2pair.tsv" # Input file
input_file="${datadir}/annotation/${POP}/final-${POP}-rv.mergedannotation.plusN2pair.Pvalthresbased.tsv" # Input file
output_prefix="model_"$model"_number_of_dimensions_"$number_of_dimensions

# RUN WATERSHED EVALUATION
echo -e "\n\n**** WATERSHED EVALUATION STARTS ****"
echo -e "\nWatershed evaluation of: $input_file" 
#Rscript evaluate_watershed.R --input $input_file --number_dimensions $number_of_dimensions --output_prefix $output_prefix --model_name $model --n2_pair_pvalue_fraction 0.009
Rscript $codedir/evaluate_watershed.R --input $input_file --number_dimensions $number_of_dimensions --output_prefix $output_prefix --model_name $model --n2_pair_pvalue_fraction 0.002
#Rscript $codedir/evaluate_watershed.R --input $input_file --number_dimensions $number_of_dimensions --output_prefix $output_prefix --model_name $model


#######################
# Run 'predict_watershed.R'
# Note for convenience, the training file is the same as the prediction file. This does not necessarily have to be the case
#######################

# Run using Watershed approximate inference
# If the number of dimensions (E) is less than equal to 4, exact inference should be used
model="RIVER"  # Can take on "Watershed_exact", "Watershed_approximate", "RIVER"
model="Watershed_approximate"  # Can take on "Watershed_exact", "Watershed_approximate", "RIVER"
model="Watershed_exact"  # Can take on "Watershed_exact", "Watershed_approximate", "RIVER"
number_of_dimensions="3" # Can take on any real number greater than or equal to one
number_of_dimensions="1" # Can take on any real number greater than or equal to one


datadir="/scratch16/abattle4/surya/datasets/WatershedAFR/data"
input_file="${datadir}/annotation/${POP}/final-${POP}-rv.mergedannotation.plusN2pair.tsv" # Input file
input_file="${datadir}/annotation/${POP}/final-${POP}-rv.mergedannotation.plusN2pair.Pvalthresbased.tsv" # Input file
output_prefix="model_"$model"_number_of_dimensions_"$number_of_dimensions

# RUN WATERSHED PREDICTION
echo -e "\n\n\n**** WATERSHED PREDICTION STARTS ****"
echo -e "\nWatershed prediction of: $input_file" 
#Rscript $codedir/predict_watershed.R --training_input $input_file --prediction_input $input_file --number_dimensions $number_of_dimensions --output_prefix $output_prefix --model_name $model
Rscript $codedir/predict_watershed.R --training_input $input_file --prediction_input $input_file --number_dimensions $number_of_dimensions --output_prefix $output_prefix --model_name $model

