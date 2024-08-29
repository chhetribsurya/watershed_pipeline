#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { watershed_modelrun } from "../modules/watershed_modelrun" params(params)

workflow WATERSHED_MODELRUN {
    take:
    pops_ch
    merged_annotfile_ch

    main:
    if (params.analysis == "GLOBAL") {
        
        // Combine the final outputs into a single channel
        combined_watershed_ch = pops_ch
            .merge(merged_annotfile_ch)
            .map{ 
                pops, merged_annotfile -> 
                tuple(pops, merged_annotfile)
        }


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

