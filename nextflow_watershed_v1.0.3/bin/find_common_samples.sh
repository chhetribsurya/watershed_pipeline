#!/bin/bash

# Function to find common samples
find_common_samples() {
    local tpm_file=$1
    local readcount_file=$2
    local covariate_file=$3
    local vcf_file=$4
    local subjids_file=$5
    local output_dir=$6
    #local pop=$7

    # Prepare the output directory
    mkdir -p "$output_dir"

    # Extract and save sample IDs
    echo -e "\nExtracting sample IDs from VCF..."
    local vcf_samples=$(bcftools query -l "$vcf_file")
    echo "Found $(echo "$vcf_samples" | wc -l) samples in VCF."
    printf "%s\n" "${vcf_samples[@]}" > "$output_dir/vcf_samples.txt"
    echo ""

    echo "Reading sample IDs from TPM file..."
    local tpm_samples=$(head -n 1 "$tpm_file" | tr '\t' '\n' | tail -n +2)
    echo "Found $(echo "$tpm_samples" | wc -l) samples in TPM file."
    printf "%s\n" "${tpm_samples[@]}" > "$output_dir/tpm_samples.txt"
    echo ""

    echo "Reading sample IDs from Read Count file..."
    local readcount_samples=$(head -n 1 "$readcount_file" | tr '\t' '\n' | tail -n +2)
    echo "Found $(echo "$readcount_samples" | wc -l) samples in Read Count file."
    printf "%s\n" "${readcount_samples[@]}" > "$output_dir/readcount_samples.txt"
    echo ""

    echo "Reading sample IDs from Covariate file..."
    local covariate_samples=$(head -n 1 "$covariate_file" | tr '\t' '\n' | tail -n +2)
    echo "Found $(echo "$covariate_samples" | wc -l) samples in Covariate file."
    printf "%s\n" "${covariate_samples[@]}" > "$output_dir/covariate_samples.txt"
    echo ""

    echo "Reading sample IDs from Subject ID file..."
    local subjid_samples=$(cut -f1 "$subjids_file")
    echo "Found $(echo "$subjid_samples" | wc -l) samples in Subject ID file."
    printf "%s\n" "${subjid_samples[@]}" > "$output_dir/subjid_samples.txt"
    echo ""

    # Finding common samples
    echo "Finding common samples..."
    comm -12 <(printf "%s\n" "${tpm_samples[@]}" | sort) <(printf "%s\n" "${readcount_samples[@]}" | sort) > step1.txt
    comm -12 <(printf "%s\n" "${covariate_samples[@]}" | sort) step1.txt > step2.txt
    comm -12 <(printf "%s\n" "${vcf_samples[@]}" | sort) step2.txt > step3.txt
    comm -12 <(printf "%s\n" "${subjid_samples[@]}" | sort) step3.txt > "$output_dir/common_samples.txt"

    local num_common=$(cat "$output_dir/common_samples.txt" | wc -l)

    if [ "$num_common" -eq 0 ]; then
        echo -e "No common samples found, exiting the pipeline run...\n\n"\
                "Please check your files for overlapping samples to ensure there are common samples"\
                "to analyze on with the pipeline.\n\n"
        exit 1
    else
        echo -e "Found $num_common common samples for analysis.\n\n"\
                "The pipeline will continue with $num_common samples. Please check the output directory"\
                " $output_dir for all relevant sample IDs including the common sample IDs.\n\n"
    fi
}

# Script execution
# Pass the file paths and parameters to the function
find_common_samples "$1" "$2" "$3" "$4" "$5" "$6" "$7"

