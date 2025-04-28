#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { generate_enrichment_plots; run_watershed_aupr } from "../modules/generate_plots" params(params)

workflow CALL_ENRICH_PLOTS {
    take:
    pops_ch
    rvfile_ch
    pval_exp_outlier_file_channel
    zscore_exp_outlier_file_channel
    
    main:
    if (params.analysis == "GLOBAL") {
        
        // Combine the final outputs into a single channel
        combined_enrichment_ch = pops_ch
            .merge(rvfile_ch, pval_exp_outlier_file_channel, zscore_exp_outlier_file_channel)
            .map{ 
                pops, rvfile, pval_expoutlier_file, zscore_expoutlier_file -> 
                tuple(pops, rvfile, pval_expoutlier_file, zscore_expoutlier_file)
        }


        // Generate enrichment plot
        generate_enrichment_plots(combined_enrichment_ch)

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
    enrichment_results = generate_enrichment_plots.out.enrichment_plots
    //auc_results = generate_enrichment_plots.out.auc_plots
}


workflow CALL_AUPR_PLOTS {
    take:
    pops_ch
    evalrds_ch
    predrds_ch
    annotfile_ch

    main:
    //def aupr_plots_enabled = (params.aupr_plots == 'true')

    // Check if AUPR plot generation is enabled
    if (params.aupr_plots) {
        log.info "AUPR plot generation enabled: ${params.aupr_plots}"

        // Combine the final outputs into a single channel
        combined_aupr_ch = pops_ch
            .merge(evalrds_ch, predrds_ch, annotfile_ch)
            .map {
                pops, evalrds_file, predrds_file, annotfile ->
                tuple(pops, evalrds_file, predrds_file, annotfile)
            }

        // Run the AUPR plot generation process
        run_watershed_aupr(combined_aupr_ch)
    
        // Emit the AUPR results
        emit:
        aupr_results = run_watershed_aupr.out.aupr_plots

    } else {
        // If AUPR plot generation is disabled, log a quiet message
        log.info "AUPR plot generation is disabled (params.aupr_plots = false). Skipping plot generation."
    }

}
