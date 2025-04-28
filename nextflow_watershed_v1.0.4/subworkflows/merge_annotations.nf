#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { merge_annotations_and_sortN2pair } from "../modules/merge_annotations" params(params)

workflow MERGE_ANNOTATIONS {
    take:
    pops_ch
    cadd_annotfile_ch
    gencode_annotfile_ch
    afreq_annotfile_ch
    ucsc_phylop_annotfile_ch
    vep_annotfile_ch
    outlier_refbasefile_ch
    n2pair_file_channel


    main:
    if (params.analysis == "GLOBAL") {
        
        // Combine the final outputs into a single channel
        combined_mergeannot_ch = pops_ch
            .merge(cadd_annotfile_ch, gencode_annotfile_ch, afreq_annotfile_ch, ucsc_phylop_annotfile_ch, vep_annotfile_ch, outlier_refbasefile_ch, n2pair_file_channel)
            .map{ 
                pops, cadd_annotfile, gencode_annotfile, afreq_annotfile, ucsc_phylop_annotfile, vep_annotfile, outlier_refbasefile, n2pair_file -> 
                tuple(pops, cadd_annotfile, gencode_annotfile, afreq_annotfile, ucsc_phylop_annotfile, vep_annotfile, outlier_refbasefile, n2pair_file)
        }

        // Merge annotation and sort N2 pairs
        merge_annotations_and_sortN2pair(combined_mergeannot_ch)



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
    merged_unscaledannot_file = merge_annotations_and_sortN2pair.out.merged_annotation_plusN2pair
    merged_scaledannot_file = merge_annotations_and_sortN2pair.out.merged_stdscaled_annotation_plusN2pair
    //n2pair_info_files = merge_annotations_and_sortN2pair.out.all
    //afreq_collapsed_file = collapse_afreq_annotation.out.collapse_out
    //afreq_uncollapsed_file = collapse_afreq_annotation.out.uncollapse_out
}

