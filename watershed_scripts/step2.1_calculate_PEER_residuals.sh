#!/bin/bash

set -o nounset -o errexit -o pipefail

## Calculate PEER factors for each tissue.
## The nubmer of PEER factors is determined by the number of samples in the tissue.
## 15 factors for < 150 samples; 30 factors for between 150 and 250 samples; 35 factors for > 250 samples

#while getopts p:e:c:g:t: flag
while getopts p:c: flag
do
  case "${flag}" in
    p) peerdir=${OPTARG};;
    c) covariates=${OPTARG};;
  esac
done

# datadir=/scratch/groups/abattle4/victor/WatershedAFR/data
# peerdir=${datadir}/data_prep/PEER
# gtex_v8_eqtl_dir=/scratch/groups/abattle4/victor/WatershedAFR/raw_data/GTEx/GTEx_Analysis_v8_eQTL
# covariates=${datadir}/data_prep/gtex_v8_eQTL_covariates.txt
# eQTL_geno=${datadir}/data_prep/gtex_2017-06-05_v8_genotypes_cis_eQTLs_012_processed.txt
# dockerimage=/home-net/home-4/xwang145@jhu.edu/code/PEER/peer-1.3.simg

#scriptdir=`dirname \$(readlink -f "\$0")`

# Get the conda hook in this shell
#eval "$(conda shell.bash hook)"
#conda activate r-env
ml  r/4.2.0
# Source conda
#source /path/to/your/conda/etc/profile.d/conda.sh

# Activate your conda environment
#conda activate r-env

runPeer() {
    traitsFileName=$1
    prefix=${traitsFileName%.log2.ztrans.txt}
    tissue=`basename $prefix`
    nsamples=$(cat $traitsFileName | wc -l) # this is actually n samples + 1
    if [ $nsamples -le 150 ]; then
        maxFactorsN=15
    elif [ $nsamples -le 249 ]; then
        maxFactorsN=30
    elif [ $nsamples -le 349 ]; then
        maxFactorsN=45
    else
        maxFactorsN=60
    fi
    outdir=${prefix}_Factors"$maxFactorsN"
    indir=${prefix}_Factors"$maxFactorsN"

    OUTFILE=${prefix}.peer.v3ciseQTLs.ztrans.txt
    if [ -f "$OUTFILE" ]; then
      echo "$OUTFILE exists."
      return
    else 
      # computing residuals
      echo "Here's the script"
      scriptdir="/scratch16/abattle4/surya/datasets/WatershedAFR/WatershedAFR/code/preprocessing/data_prep" 
      Rscript ${scriptdir}/step2.1_calculate_PEER_residuals.R $traitsFileName $covariates ${indir}/factors.tsv ${OUTFILE} &> ${outdir}/log.residuals.txt

      #Rscript ${scriptdir}/step2.1_calculate_PEER_residuals.R $traitsFileName $covariates \
      #  ${indir}/factors.tsv ${gtex_v8_eqtl_dir}/${tissue}.v8.egenes.txt.gz \
      #  $eQTL_geno ${prefix}.peer.v3ciseQTLs.ztrans.txt &> ${outdir}/log.residuals.txt  
    fi

}

export -f runPeer
export peerdir
#export gtex_v8_eqtl_dir
export covariates
#export eQTL_geno
#export TMPDIR
#export scriptdir

#parallel --jobs 3 runPeer ::: ${peerdir}/*.log2.ztrans.txt
for exprfile in ${peerdir}/*.log2.ztrans.txt;do
    echo -e "\nProcessing expression file: $exprfile ...\n";
    runPeer $exprfile
done

echo "PEER FACTOR RESIDUAL COMPUTE DONE!"
