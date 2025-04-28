#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process identify_samples_to_analyze {
    debug true
    tag { "JOB ${pop}" }
    publishDir "$params.out_dir/outputs/", mode : 'copy'

    input:
        tuple val(subjids_file), val(pop)

    output:
       path("sample_select/$pop/common_samples.txt"), emit: samples
       path("sample_select/$pop/*.txt"), emit: all

    script:
        def outputdir = "sample_select/$pop"

        """
        mkdir -p $outputdir
        echo -e "\n\n"
        echo "---------------------------------------------------"
        bash ${launchDir}/bin/find_common_samples.sh $params.tpm_infile $params.readcount_infile $params.covariate_infile $params.rv_file "$subjids_file" "$outputdir"
        echo ""
        echo "---------------------------------------------------"
        echo -e "\n\n"
        """
}
