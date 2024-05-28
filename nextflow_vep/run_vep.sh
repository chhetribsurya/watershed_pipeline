#!/bin/bash

# Nextflow sample run with path defined 
# NXF_PATH="/scratch16/abattle4/surya/tools/nextflow_run/"
# ${NXF_PATH}/nextflow run vep_finaltest.nf -c vep_finaltest.config

# Nextflow Sample run: example usage 
# nextflow main_chromwise.nf --skip_cache --resultoutdir /scratch16/abattle4/surya/datasets/WatershedAFR/data/annotation --inputvcf gene-GLOBAL-rv.CADD.vcf -profile local

#*****************************************************************************************************************************************************************************************

# Nextflow VEP Annotation Pipeline

#*****************************************************************************************************************************************************************************************

# Population 
POP=$1
POP="EUR"

# Load java module
ml java/19

# Set up directory
outdir="/scratch16/abattle4/surya/datasets/WatershedAFR/data/annotation/${POP}"
chromwise_outdir="/scratch16/abattle4/surya/datasets/WatershedAFR/data/annotation/${POP}"
vcf_file="/scratch16/abattle4/surya/datasets/WatershedAFR/data/rare_variants_gnomad/${POP}/gene-${POP}-rv.CADD.vcf"
vcfchromfileDir="/scratch16/abattle4/surya/datasets/WatershedAFR/data/rare_variants_gnomad/${POP}/chromwise"

# Define Nextflow scripts containing dir
nxf_dir="/scratch16/abattle4/surya/datasets/watershed_rarevars/project_final"

# FOR CHROMWISE VEP RUN
# Skips cache buildup
nextflow ${nxf_dir}/main_chromwise_final.nf -c ${nxf_dir}/nextflow.config --skip_cache --resultoutdir $chromwise_outdir --vcfDir $vcfchromfileDir -profile slurm

# Includes Cache buidup
#nextflow ${nxf_dir}/main_chromwise_final.nf -c ${nxf_dir}/nextflow.config --resultoutdir $chromwise_outdir --vcfDir $vcfchromfileDir -profile slurm

# FOR WHOLE-CHROMOSOME BASED RUN
#nextflow ${nxf_dir}/main_nonchromwise.nf --skip_cache --resultoutdir $outdir --inputvcf $vcf_file -profile slurm

