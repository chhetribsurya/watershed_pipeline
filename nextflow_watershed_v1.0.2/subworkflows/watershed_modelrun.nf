#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*include { process_data; call_phylopscore_annotation } from "../modules/call_annotation_process.nf" params(params)*/

process watershed_modelrun {
    debug true
    tag { "JOB ${pop}" }
    //label 'watershed_py_env'
    publishDir "$params.out_dir/outputs/7_watershed_results/", mode : 'copy'

    input:
        tuple val(pop), val(merged_annotation_file)

    output:
        path "Watershed/model_outputs/${pop}/${pop}*evaluation_object.rds", emit: evaluation
        path "Watershed/model_outputs/${pop}/${pop}*posterior_probability.txt", emit: posterior
        path "Watershed/model_outputs/${pop}/${pop}*prediction_object.rds", emit: prediction
        path "logs/${pop}/watershed.out", optional: true, emit: outfile

    script:
        def model = "Watershed_exact" // [RIVER, Watershed_approximate]
        def num_of_dim = 1
        def outputdir = "model_outputs/${pop}" 
        def outprefix = "${pop}_model_${model}_dim_${num_of_dim}" //output file
        def watershed_dir = "${launchDir}/src/Watershed"
        def logdir = "logs/${pop}"

        """
        #!/usr/bin/bash
        mkdir -p $logdir
        cp -r $watershed_dir .
        # Set the correct working directory for the R script
        cd Watershed
        mkdir -p $outputdir
        echo "Here's the content in the Watershed directory:"
        ls -la
        echo -e "\n***** Watershed model prediction for POP: ${pop}, MERGED Annotation file: $merged_annotation_file, MODEL: $model, NUM OF DIM: $num_of_dim, OUTPUT DIR: $outputdir, OUTPUT FILE PREFIX: $outprefix, WATERSHED SRC DIR: $watershed_dir"
        bash "${launchDir}/bin/step17_watershed_model_call-final.sh" \
            "$pop" \
            "${merged_annotation_file}" \
            "$outputdir" \
            "$outprefix" \
            "$num_of_dim" \
            "$model" \
            "Watershed"
        echo -e "\n\n**** PIPELINE SUCCESSFULLY EXECUTED! ****"
        #ls -ltr
        """
}


workflow WATERSHED_MODELRUN {
    take:
    pops_ch
    merged_annotfile_ch

    main:
    if (params.analysis == "GLOBAL") {
        
        //Channel.value("$params.cohort_name").set { pops_ch }
        //pops_ch.merge(rv_file_ch, vep_loftee_file_ch)

        // Combine the final outputs into a single channel
        combined_watershed_ch = pops_ch
            .merge(merged_annotfile_ch)
            .map{ 
                pops, merged_annotfile -> 
                tuple(pops, merged_annotfile)
        }

        combined_watershed_ch.view({ "SUBWORKFLOW WATERSHED CHANNEL combined output: $it" })

        // Run watershed model prediction
        watershed_modelrun(combined_watershed_ch)



       /*------------------------------------------------------------------*/

    } else if (params.analysis == "STRATIFIED") {

        // Joining channels based on the key
        combined_ch = residual_expr_keyed
            .join(gene_rarevar_keyed)
            .map { key, residual_expr, gene_rarevar ->
                tuple(key, residual_expr, gene_rarevar)
            }
            .filter { key, residual_expr, gene_rarevar -> key == "EUR" || key == "AFR" }
            .view({ key, residual_expr, gene_rarevar -> "Final generare-variant combined tuple: ($key, $residual_expr, $gene_rarevar)" })
        
        // Call expression outliers using the synchronized combined channel
        call_expr_outliers(combined_ch)
        call_expr_outliers.out.outliers.view({"Here's the outlier file list: $it"})
        call_expr_outliers.out.pops_outliers.view({"Here's the outlier pop: $it"})

    }

    emit:
    watershed_evaluation = watershed_modelrun.out.evaluation
    watershed_posterior = watershed_modelrun.out.posterior
    watershed_prediction = watershed_modelrun.out.prediction
    //n2pair_info_files = merge_annotations_and_sortN2pair.out.all
    //afreq_collapsed_file = collapse_afreq_annotation.out.collapse_out
    //afreq_uncollapsed_file = collapse_afreq_annotation.out.uncollapse_out
}

