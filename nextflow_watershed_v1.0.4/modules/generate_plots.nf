#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process generate_enrichment_plots {
    tag { "JOB ${pop}" }
    publishDir "$params.out_dir/outputs/8_plots/", mode: 'copy'

    input:
        tuple val(pop), val(rv_file), val(pval_exp_outlier_file), val(zscore_exp_outlier_file)

    output:
        //path "enrichment_plots/${pop}/*.pdf", emit: enrichment_plots
        path "enrichment_plots/${pop}/rarevar_enrichment_${pop}/*.pdf", emit: enrichment_plots

    script:
        def outputdir = "enrichment_plots/${pop}"
        
        /* Define R scripts based on the methods */
        def r_script_pval = "${launchDir}/bin/step7_rare_var_enrichment-pval.R"
        def r_script_zscore = "${launchDir}/bin/step7_rare_var_enrichment.R"

        /* Handle enrichment method */
        if (params.enrichment_method == "pvalue") {
            """
            mkdir -p $outputdir
            echo -e "\nEnrichment method used: pvalue"
            echo -e "Script used for enrichment run: $r_script_pval\n"
            Rscript "$r_script_pval" \
                "$outputdir" \
                "$pop" \
                "${pval_exp_outlier_file}" \
                "${rv_file}"
            echo -e "\n\n*** Enrichment plot generated for pval version ***"
            """
        } else if (params.enrichment_method == "zscore") {
            """
            mkdir -p $outputdir
            echo -e "\nEnrichment method used: zscore"
            echo -e "Script used for enrichment run: $r_script_zscore\n"
            Rscript "$r_script_zscore" \
                "$outputdir" \
                "$pop" \
                "${zscore_exp_outlier_file}" \
                "${rv_file}"
            echo -e "\n\n*** Enrichment plot generated for zscore version***"
            """
        } else if (params.enrichment_method == "both") {
            """
            mkdir -p $outputdir
            echo -e "\nEnrichment method used: both"
            
            # Run the pvalue method
            echo -e "Running pvalue enrichment script: $r_script_pval\n"
            Rscript "$r_script_pval" \
                "$outputdir" \
                "$pop" \
                "${pval_exp_outlier_file}" \
                "${rv_file}"
            
            # Run the zscore method
            echo -e "Running zscore enrichment script: $r_script_zscore\n"
            Rscript "$r_script_zscore" \
                "$outputdir" \
                "$pop" \
                "${zscore_exp_outlier_file}" \
                "${rv_file}"
            
            echo -e "\n\n*** Enrichment plots generated for both methods ***"
            """
        } else {
            error "Invalid enrichment method specified: ${params.enrichment_method}. Valid options are 'pvalue', 'zscore', or 'both'."
        }
}

process run_watershed_aupr {
    tag { "JOB ${pop}" }
    publishDir "$params.out_dir/outputs/8_plots/", mode: 'copy'

    input:
        tuple val(pop), path(eval_rds), path(pred_rds), path(annot_file)

    output:
        path "aupr_plots/${pop}/*.pdf", emit: aupr_plots
        //path "aupr_plots/${pop}/*.csv", emit: aupr_files

    script:
        def output_dir = "aupr_plots/${pop}"

        """
        echo -e "Running aupr analysis\n"
        mkdir -p $output_dir
        Rscript "${launchDir}/bin/step18_watershed_aupr.R" \
            "${eval_rds}" \
            "${pred_rds}" \
            "${annot_file}" \
            "${output_dir}" \
            "${pop}"
        echo -e "\n\n*** AUPR plots generated for watershed results ***"
        """
}

