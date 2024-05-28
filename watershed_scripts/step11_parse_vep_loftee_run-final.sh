#!/bin/bash

#SBATCH --job-name=ParseVEPannotation
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=20G

# Population variable from the argument
POP="$1"

# Function to display usage instructions
function usage {
    echo "Usage: $0 <POP variable>"
    echo "Performs the following steps:"
    echo "1) Activates a conda environment."
    echo "2) Processes chromosome-wise VCF files (e.g., gene-EUR-rv.chr*.vep.loftee.vcf)."
    echo "3) bgzip's each VCF file while retaining the original."
    echo "4) Indexes each gzipped VCF with tabix."
    echo "5) Runs a Python script using the gzipped VCF files as input."
    echo "6) Consolidates the output into the specified TSV file."
    exit 1
}

# Check if the correct number of arguments is provided
#if [ $# -ne 2 ]; then
#    usage
#fi

# Check if user asks for help
if [[ "$1" == "--help" || "$1" == "-h" ]]; then
    usage
fi

# Step 1: Activate conda environment
watershed="/home/schhetr1/anaconda3/envs/watershed"
eval "$(conda shell.bash hook)"
conda activate $watershed

if [[ $? -ne 0 ]]; then
    echo "Failed to activate watershed conda environment. Exiting."
    exit 1
fi

# File paths
datadir1="/scratch16/abattle4/surya/datasets/WatershedAFR/data"
codedir1="/scratch16/abattle4/surya/datasets/WatershedAFR/WatershedAFR/code/preprocessing/annotation"
outputfile="${datadir1}/annotation/${POP}/vep/gene-${POP}-rv.vep.loftee.snv.tsv"
vep_loftee_prefix="${datadir1}/annotation/${POP}/vep/gene-${POP}-rv"

# Check if output file already exists
if [ -f "$outputfile" ]; then
    echo "Output file already exists: $outputfile"
    echo "Skipping Python script execution."
    exit 0
fi

# Loop through all chromosome-specific VCF files and process them
for file in ${vep_loftee_prefix}.chr*.vep.loftee.vcf; do
    if [ ! -f "$file" ]; then
        echo "File not found: $file"
        continue
    fi
 
    gz_file="${file}.gz"
    if [ -f "$gz_file" ]; then
        echo "Bgzipped file already exists: $gz_file"
    else
        echo -e "\nBgzipping file: $file"
        bgzip -c "$file" > "$gz_file"
        echo -e "\nBgzip Done ..."
    fi

    if [ -f "${gz_file}.tbi" ]; then
        echo "Tabix index already exists for: $gz_file"
    else
        echo -e "\nTabixing file: $gz_file"
        tabix -p vcf "$gz_file"
        echo -e "\nTabix Done ..."
    fi

done


# Run the Python script with the annotation prefix and output file
echo -e "\nPython parsing of VEP VCF files ..."
python ${codedir1}/step11.1_parse_vep_loftee_chromwise-final.py --anno_prefix $vep_loftee_prefix --outputfile $outputfile
#python ${codedir1}/step11.1_parse_vep_loftee_chromwise.py --anno_prefix $vep_loftee_prefix --outputfile $outputfile
#python ${codedir1}/step11.1_parse_vep_loftee_Fullchrom.py --anno_prefix $vep_loftee_prefix --outputfile $outputfile
echo -e "\n\n*** Steps completed successfully *** \n\n"
