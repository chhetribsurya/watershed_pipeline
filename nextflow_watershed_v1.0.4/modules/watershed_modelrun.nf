#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process watershed_modelrun {
    tag { "JOB ${pop}" }
    publishDir "$params.out_dir/outputs/7_watershed_results/", mode : 'copy'

    input:
        tuple val(pop), val(merged_annotation_file), val(pval_thresh_file)

    output:
        path "Watershed/model_outputs/${pop}/${pop}*evaluation_object.rds", emit: evaluation
        path "Watershed/model_outputs/${pop}/${pop}*posterior_probability.txt", emit: posterior
        path "Watershed/model_outputs/${pop}/${pop}*prediction_object.rds", emit: prediction
        path "logs/${pop}/watershed.log", emit: outfile

    script:
        def model = "Watershed_exact" // [RIVER, Watershed_approximate]
        def num_of_dim = 1
        def outputdir = "model_outputs/${pop}" 
        def outprefix = "${pop}_model_${model}_dim_${num_of_dim}" //output file
        def watershed_dir = "${launchDir}/src/Watershed"
        def logdir = "logs/${pop}"

        // Conditionally select the bash script based on params.outlier_method
        def bash_script = params.outlier_method == "zscore" 
            ? "${launchDir}/bin/step17_watershed_model_call-final-zscoreversion.sh"
            : "${launchDir}/bin/step17_watershed_model_call-final-test.sh"

        """
        mkdir -p $logdir
        cp -r $watershed_dir .
        cd Watershed
        mkdir -p $outputdir
        echo -e "\nOutlier method used: $params.outlier_method"
        echo -e "Script used for watershed run: $bash_script\n"
        bash "$bash_script" \
            "$pop" \
            "${merged_annotation_file}" \
            "$outputdir" \
            "$outprefix" \
            "$num_of_dim" \
            "$model" \
            "Watershed" \
            "${pval_thresh_file}"
        cd ..
        cp \$PWD/.command.out "$logdir/watershed.log"
        echo -e "\n\n**** PIPELINE SUCCESSFULLY EXECUTED! ****"
        """
}
