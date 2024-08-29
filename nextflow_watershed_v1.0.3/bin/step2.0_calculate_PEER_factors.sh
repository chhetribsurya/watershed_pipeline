#!/bin/bash

set -o nounset -o errexit -o pipefail

## Calculate PEER factors for each tissue.
## The nubmer of PEER factors is determined by the number of samples in the tissue.
## 15 factors for < 150 samples; 30 factors for between 150 and 250 samples; 35 factors for > 250 samples


#while getopts p:c:t:d: flag
#while getopts p:s:t:d: flag
#while getopts p:s:d: flag
while getopts p:s: flag
do
  case "${flag}" in
    p) peer_infile=${OPTARG};;
    s) scriptdir=${OPTARG};;
  esac
done

#[ ! -d "$tmpdir" ] && mkdir -p "$tmpdir"
#export TMPDIR=$tmpdir

runPeer() {
    traitsFileName=$1
    scriptdir=$2
    #dockerimage=$3
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
    maxIterations=10000
    boundTol=0.001
    varTol=0.00001
    e_pa=0.1
    e_pb=10
    a_pa=0.001
    a_pb=0.1
    outdir=${prefix}_Factors"$maxFactorsN"
    indir=${prefix}_Factors"$maxFactorsN"
    echo $outdir

    mkdir -p $outdir

    ## actual calculation of peer factors
    echo "Rscript step2.0_calculate_PEER_factors.R $traitsFileName $maxFactorsN $maxIterations $boundTol $varTol $e_pa $e_pb $a_pa $a_pb $outdir $tissue" > ${outdir}/log.txt

    # Bind all directories used
    #bindpaths="$scriptdir,$(dirname $traitsFileName)"
 
    Rscript ${scriptdir}/step2.0_calculate_PEER_factors.R $traitsFileName $maxFactorsN $maxIterations $boundTol $varTol $e_pa $e_pb $a_pa $a_pb $outdir $tissue >> ${outdir}/log.txt 2>&1 
    
}

export -f runPeer

# INPUT Parameters
export peer_infile
export scriptdir
#export dockerimage
#export TMPDIR

#parallel --jobs 3 runPeer ::: ${peerdir}/*.log2.ztrans.txt
#for exprfile in ${peerdir}/*.log2.ztrans.txt;do
echo -e "\nProcessing expression file: $peer_infile ...\n";
#runPeer $peer_infile $scriptdir $dockerimage
runPeer $peer_infile $scriptdir
#done

echo "PEER FACTOR COMPUTE DONE!"

