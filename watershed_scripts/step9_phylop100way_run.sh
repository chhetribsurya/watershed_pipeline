datadir1="/scratch16/abattle4/surya/datasets/WatershedAFR/data"
codedir="/scratch16/abattle4/surya/datasets/WatershedAFR/WatershedAFR/code/preprocessing/annotation"
annodir1="${datadir1}/annotation/ucsc"
chromwise_dir="${datadir1}/rare_variants_gnomad/chromwise"

ucsc_path="/scratch16/abattle4/surya/tools/ucsc"
bigwig="/scratch16/abattle4/surya/tools/ucsc/phylop_scores/hg38.phyloP100way.bw"
bedgraph="${annodir1}/hg38.phyloP100way.bedGraph"

# Create ucsc annotation directory if it doesn't exist
mkdir -p $annodir1

# # RUN bigwig to bedgraph
# ${ucsc_path}/bigWigToBedGraph $bigwig $bedgraph

# Load R module
ml r/4.2.0 

## For each chromosome-wise file, convert rare variants to bed format
#for file in ${chromwise_dir}/gene-GLOBAL-rv.*.txt; do
#    echo -e "\nProcessing file: $file \n"
#    # Convert rare variants to bed format
#    Rscript $codedir1/step8.1_rare_variants_to_bed.R --RV ${file}
#    
#done

rarevar_file="${datadir1}/rare_variants_gnomad/gene-GLOBAL-rv.txt"
echo -e "\nProcessing rare variant file: $rarevar_file \n"
# Convert rare variants to bed format
#Rscript $odedir1/step8.1_rare_variants_to_bed.R --RV ${rarevar_file}
bash ${codedir}/step9.1_phylop100way.sh -b $bedgraph -r ${datadir1}/rare_variants_gnomad/gene-GLOBAL-rv.bed GLOBAL
