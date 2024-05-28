#/usr/bin/env bash

# Directory settings
datadir1="/scratch16/abattle4/surya/datasets/WatershedAFR/data"
codedir1="/scratch16/abattle4/surya/datasets/WatershedAFR/WatershedAFR/code/preprocessing/annotation"

chromwise_dir="${datadir1}/rare_variants_gnomad/chromwise"
cadd_loc="/scratch16/abattle4/surya/tools/CADD-scripts"
pop="GLOBAL"
corenum=6

log_dir="/scratch16/abattle4/surya/datasets/WatershedAFR/WatershedAFR/code/preprocessing/annotation/logfiles"
outdir="/scratch16/abattle4/surya/datasets/WatershedAFR/data/annotation/cadd"

# Create logfiles directory if it doesn't exist
mkdir -p "${log_dir}"
mkdir -p "${outdir}"

for file in ${chromwise_dir}/gene-GLOBAL-rv.*.vcf; do
    chrom=$(basename $file | cut -d'.' -f2)
    echo -e "\nProcessing chrom: $chrom"
    
    sbatch \
    --job-name=CADD_annotate_${chrom} \
    --output="${log_dir}/CADD_annotate_${chrom}.log" \
    --time=15:00:00 \
    --mem=15G \
    --cpus-per-task=${corenum} \
    --wrap="${codedir1}/step7.1_cadd_annotate.sh ${corenum} ${chrom} ${pop} ${file} ${outdir}"
done
