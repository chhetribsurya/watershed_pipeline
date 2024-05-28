#!/usr/bin/env bash

# Array of population variables
POPULATIONS=("AFR")
POPULATIONS=("EUR")
POPULATIONS=("GLOBAL" "EUR" "AFR")

# CADD Annotation
echo -e "\n\n***** CADD Annotation starts: $POP *****\n"
for POP in "${POPULATIONS[@]}"; do
    # CADD Pre-processing of rare variant to VCF chromwise file for CADD run
    echo -e "*** CADD Pre-Processing of Pop: $POP ***"
    bash step6_chromwise_rare_variants_to_vcf.sh "$POP"
    
    # RUN CADD annotation chromwise
    echo -e "*** CADD Processing of Pop: $POP ***"
    bash step7_cadd_slurm_run.sh "$POP"
     
done

# UCSC Annotation
echo -e "\n\n***** UCSC Annotation starts: $POP *****\n"
for POP in "${POPULATIONS[@]}"; do
    # UCSC pre-processing of rarevar file to bedfile for phyloP query format
    echo -e "*** UCSC-phyloP Pre-Processing of Pop: $POP ***"
    job1=$(sbatch --parsable step8_rare_variants_to_bed.sh "$POP")
    sbatch step8_rare_variants_to_bed.sh "$POP"
    
    # Run UCSC annotation for phyloP scoring
    echo -e "*** UCSC-phyloP Processing of Pop: $POP ***"
    sbatch --dependency=afterok:$job1 step9_phylop100way_run.sh $POP
    ##sbatch step9_phylop100way_run.sh $POP
done


# VEP Annotation
echo -e "\n\n***** VEP Annotation starts: $POP *****\n"
for POP in "${POPULATIONS[@]}"; do
    # VEP Pre-processing of rare variant to VCF chromwise file for VEP run
    # Exactly similar as CADD based pre-processing. Thus below run not required

    ##echo -e "\n*** VEP Pre-Processing of Pop: $POP ***"
    ##bash step6_chromwise_rare_variants_to_vcf.sh "$POP"

    # Run VEP annotation chromwise
    echo -e "*** VEP Processing of Pop: $POP ***"
    bash vep_final.sh $POP
    ##bash vep_finaltest.sh $POP

    # Parse VEP annotation
    echo -e "*** VEP file parsing of Pop: $POP ***"
    sbatch step11_parse_vep_loftee_run-final.sh $POP
done

# GENCODE based TSS TES annotation
echo -e "\n\n***** Gencode Annotation starts: $POP *****\n"
for POP in "${POPULATIONS[@]}"; do 
    # Run Gencode annotation to compute dist to TSS and TES from rarevars
    echo -e "*** Gencode Annot Processing of Pop: $POP ***"
    sbatch step10_gencode_tss_tes_distance.sh $POP
done


