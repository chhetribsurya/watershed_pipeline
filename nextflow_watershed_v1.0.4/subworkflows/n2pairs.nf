#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { generate_n2_pairs } from "../modules/n2pairs" params(params)

workflow CALL_N2_PAIRS {
    take:
    pops_ch
    rv_file_ch
    exp_outlier_file_ch

    main:
    if (params.analysis == "GLOBAL") {
        
        // Combine the final outputs into a single channel
        combined_n2pair_ch = pops_ch
            .merge(rv_file_ch, exp_outlier_file_ch)
            .map{ pops, rv_file, exp_outlier -> tuple(pops, rv_file, exp_outlier)
        }

        // Process and analyze N2pair annotation
        generate_n2_pairs(combined_n2pair_ch)



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
    n2pair_file = generate_n2_pairs.out.n2pairs
    n2pair_info_files = generate_n2_pairs.out.all
    //afreq_collapsed_file = collapse_afreq_annotation.out.collapse_out
    //afreq_uncollapsed_file = collapse_afreq_annotation.out.uncollapse_out
}

