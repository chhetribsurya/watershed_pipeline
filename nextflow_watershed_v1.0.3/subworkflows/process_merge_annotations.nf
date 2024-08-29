#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { process_merged_annotations } from "../modules/process_merge_annotations" params(params)

workflow PROCESS_MERGED_ANNOTATION {
    take:
    pops_ch
    merged_scaled_annotfile_channel

    main:
    if (params.analysis == "GLOBAL") {
        
        // Combine the final outputs into a single channel
        combined_process_mergeannot_ch = pops_ch
            .merge(merged_scaled_annotfile_channel)
            .map{ 
                pops, merged_annotfile -> 
                tuple(pops, merged_annotfile)
        }

        // Merge annotation and sort N2 pairs
       process_merged_annotations(combined_process_mergeannot_ch)



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
    processed_mergedannot_file = process_merged_annotations.out.processed_file
    featurecount_file = process_merged_annotations.out.feature_counts
    //test_channel  = combined_process_mergeannot_ch
    //n2pair_info_files = merge_annotations_and_sortN2pair.out.all
    //afreq_collapsed_file = collapse_afreq_annotation.out.collapse_out
    //afreq_uncollapsed_file = collapse_afreq_annotation.out.uncollapse_out
}

