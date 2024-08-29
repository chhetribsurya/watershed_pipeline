#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process watershed_modelrun {
    tag { "JOB ${pop}" }
    publishDir "$params.out_dir/outputs/7_watershed_results/", mode : 'copy'

    input:
        tuple val(pop), val(merged_annotation_file)

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

        """
        mkdir -p $logdir
        cp -r $watershed_dir .
        cd Watershed
        mkdir -p $outputdir
        bash "${launchDir}/bin/step17_watershed_model_call-final.sh" \
            "$pop" \
            "${merged_annotation_file}" \
            "$outputdir" \
            "$outprefix" \
            "$num_of_dim" \
            "$model" \
            "Watershed"
        cd ..
        cp \$PWD/.command.out "$logdir/watershed.log"
        echo -e "\n\n**** PIPELINE SUCCESSFULLY EXECUTED! ****"
        """
}
