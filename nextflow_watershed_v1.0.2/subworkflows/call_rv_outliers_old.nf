#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { process_data; normalize_expression; compute_peerfactors; compute_expr_residuals; call_expr_outliers, call_expr_outliers_only } from "../modules/expr_outlier_process" params(params)

workflow CALL_OUTLIERS {

    main:
    if (params.analysis == "GLOBAL") {
        
        // Process all individuals together
        input_ch = tuple(params.subjids_file, params.cohort_name)
        process_data(input_ch) 

        // specific emit based output
        process_data.out.all.view({"\nHere's the list of ALL processed output file: $it \n"}) //outputs all (everything)
        
        // Collect the outputs from process_data
        tpms_ch = process_data.out.tpms
        reads_ch = process_data.out.reads
        covariates_ch = process_data.out.covariates
        pops_ch = process_data.out.pops

        // Normalize expression process
        normalize_expression(tpms_ch, reads_ch, covariates_ch, pops_ch)

        // Compute PEER factors
        lognorm_expr_ch = normalize_expression.out.lognormalized
        pops_norm_ch = normalize_expression.out.pops_norm

        compute_peerfactors(lognorm_expr_ch, pops_norm_ch)
        compute_peerfactors.out.peer_factors.view({"Here's the PEER factor file: $it"})
        compute_peerfactors.out.pops_peer.view({"Here's the PEER factor pop: $it"})

        // Compute expression residuals
        peer_factors_ch = compute_peerfactors.out.peer_factors
        pops_peer_ch = compute_peerfactors.out.pops_peer


        // Combine the final outputs into a single channel and form tuples
        combined_ch = pops_norm_ch
            .merge(lognorm_expr_ch, peer_factors_ch, covariates_ch) { pops_norm, lognorm_expr, peer_factors, covariates ->
            tuple(pops_norm, lognorm_expr, peer_factors, covariates)
        }
        combined_ch.view({ pops_norm, lognorm_expr, peer_factors, covariates -> "Final combined output: ($pops_norm, $lognorm_expr, $peer_factors, $covariates)" })
        
        compute_expr_residuals(combined_ch)
        compute_expr_residuals.out.residuals.view({"Here's the residual factor file: $it"})
        compute_expr_residuals.out.pops_residuals.view({"Here's the residual factor pop: $it"})

        // Call expression outliers
        residual_expr_ch = compute_expr_residuals.out.residuals
        pops_peer_ch = compute_expr_residuals.out.pops_residuals

        // Combine the final outputs into a single channel and form tuples
        //Channel.fromPath($params.gene_rarevar_file).set(gene_rarevar_ch)

        // Setting rare-variant file channel
        //Channel.fromPath('./data/gene-GLOBAL-rv.txt').set { gene_rarevar_ch }
        //Channel.fromPath('./data/gene*GLOBAL*rv.txt').set { gene_rarevar_ch }

        Channel.fromPath("$params.rv_file")
            | set { gene_rarevar_ch }

        pops_peer_ch.view({"Here's the outlier pops_peer_ch: $it"})
        residual_expr_ch.view({"Here's the outlier residual_expr_ch: $it"})
        gene_rarevar_ch.view({"Here's the outlier gene_rarevar_ch: $it"})
        combined_ch = pops_peer_ch
            .merge(residual_expr_ch, gene_rarevar_ch)
            .map { pops_peer, residual_expr, gene_rarevar ->
            tuple(pops_peer, residual_expr, gene_rarevar)
        }

        combined_ch.view({ pops_peer, residual_expr, gene_rarevar  -> "Final Global outlier call step combined output: ($pops_peer, $residual_expr, $gene_rarevar)" })
        
        call_expr_outliers(combined_ch)
        call_expr_outliers.out.outliers.view({"Here's the outlier file list: $it"})
        call_expr_outliers.out.pops_outliers.view({"Here's the outlier pop: $it"})


    } else if (params.analysis == "STRATIFIED") {

        // Split the population file into individual population files and process each one separately
        split_ch = split_population(params.subjids_file).flatten()

        // Create channel
        pop_ch = split_ch.map { file -> tuple(file, file.baseName.split('_')[0]) } // extract population label
        pop_ch.view({"\nHere's the tuple based pop_ch: $it \n"})

        // Process population data
        process_data(pop_ch) 

        /* Collect the outputs from process_data */
        tpms_ch = process_data.out.tpms
        reads_ch = process_data.out.reads
        covariates_ch = process_data.out.covariates
        pops_ch = process_data.out.pops

        // Normalize expression process
        normalize_expression(tpms_ch, reads_ch, covariates_ch, pops_ch)

        lognorm_expr_ch = normalize_expression.out.lognormalized
        pops_norm_ch = normalize_expression.out.pops_norm

        // Compute PEER factors
        compute_peerfactors(lognorm_expr_ch, pops_norm_ch)

        peer_factors_ch = compute_peerfactors.out.peer_factors
        pops_peer_ch = compute_peerfactors.out.pops_peer

        compute_peerfactors.out.peer_factors.view({"Here's the PEER factor file: $it"})
        compute_peerfactors.out.pops_peer.view({"Here's the PEER factor pop: $it"})


        // Keying lognorm_expr_ch using pops_norm_ch
        lognorm_expr_ch = normalize_expression.out.lognormalized
        pops_norm_ch = normalize_expression.out.pops_norm

        lognorm_expr_keyed = lognorm_expr_ch
            .merge(pops_norm_ch) { lognorm_expr, pops_norm -> tuple(pops_norm, lognorm_expr) }
        lognorm_expr_keyed.view({ key, lognorm_expr -> "Keyed lognorm_expr_ch: ($key, $lognorm_expr)" })

        // Keying covariates_ch using pops_ch
        covariates_ch = process_data.out.covariates
        pops_ch = process_data.out.pops

        covariates_keyed = covariates_ch
            .merge(pops_ch) { covariates, pops -> tuple(pops, covariates) }
        covariates_keyed.view({ key, covariates -> "Keyed covariates_ch: ($key, $covariates)" })

        // Keying peer_factors_ch using pops_peer_ch
        peer_factors_ch = compute_peerfactors.out.peer_factors
        pops_peer_ch = compute_peerfactors.out.pops_peer

        peer_factors_keyed = peer_factors_ch
            .merge(pops_peer_ch) { peer_factors, pops_peer -> tuple(pops_peer, peer_factors) }
        peer_factors_keyed.view({ key, peer_factors -> "Keyed peer_factors_ch: ($key, $peer_factors)" })

        // Joining channels based on the key
        combined_ch = lognorm_expr_keyed
            .join(peer_factors_keyed)
            .join(covariates_keyed)
            .map { key, lognorm_expr, peer_factors, covariates ->
                tuple(key, lognorm_expr, peer_factors, covariates)
            }
            .filter { key, lognorm_expr, peer_factors, covariates -> key == "EUR" || key == "AFR" }
            .view({ key, lognorm_expr, peer_factors, covariates -> "Final combined tuple: ($key, $lognorm_expr, $peer_factors, $covariates)" })


        // Compute residuals using the synchronized combined channel
        compute_expr_residuals(combined_ch)
        compute_expr_residuals.out.residuals.view({"Here's the residual factor file: $it"})
        compute_expr_residuals.out.pops_residuals.view({"Here's the residual factor pop: $it"})
        
        // Call expression outliers
        residual_expr_ch = compute_expr_residuals.out.residuals
        pops_expr_ch = compute_expr_residuals.out.pops_residuals
        
        // Keying residual_expr_ch using pops_expr_ch
        residual_expr_keyed = residual_expr_ch
            .merge(pops_expr_ch) { residual_expr, pops_expr -> tuple(pops_expr, residual_expr) }
        residual_expr_keyed.view({ key, residual_expr -> "Keyed residual_expr_ch: ($key, $residual_expr)" })

        // Combine the final outputs into a single channel and form tuples
        // Channel.fromPath($params.gene_rarevar_file).set(gene_rarevar_ch)

        // Keying gene rare variant files with pop-names
        Channel.fromPath('./data/gene*rv.txt')
            | map { tuple( it.baseName.split('-')[1], it ) } // extract population label
            | set { gene_rarevar_keyed }

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
    outlier_pval_file  = call_expr_outliers.out.outlier_pvalthresh
    global_outlier_removed_file  = call_expr_outliers.out.outlier_removed_pval
    residual_exp_file  = compute_expr_residuals.out.residuals
    combined_ch_list = combined_ch

}
