#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*include { process_data; call_phylopscore_annotation } from "../modules/call_annotation_process.nf" params(params)*/

process process_merged_annotations {
    debug true
    tag { "JOB ${pop}" }
    //label 'watershed_py_env'
    publishDir "$params.out_dir/outputs/6_watershed_annotation/", mode : 'copy'

    input:
        tuple val(pop), val(merged_annotation_file)

    output:
        path "processed_merged_annotation/${pop}/final-${pop}-rv.mergedannotation.plusN2pair.Pvalthresbased.stdscaled.annot.clean.tsv", emit: processed_file
        //path "processed_merged_annotation/${pop}/final-${pop}-rv.mergedannotation.plusN2pair.Pvalthresbased.annot.clean.tsv", emit: processed_file
        path "processed_merged_annotation/${pop}/feature_counts_${pop}.tsv", emit: feature_counts
        path "processed_merged_annotation/${pop}/dropped_columns_${pop}.txt", emit: dropped_columns

    script:
        def annotation_dir = "processed_merged_annotation"
        def outputdir = "$annotation_dir/${pop}"
        def output_file = "$outputdir/final-${pop}-rv.mergedannotation.plusN2pair.Pvalthresbased.stdscaled.annot.clean.tsv"
        //def output_file = "$outputdir/final-${pop}-rv.mergedannotation.plusN2pair.Pvalthresbased.clean.tsv"

        """
        #!/usr/bin/bash
        mkdir -p $outputdir
        echo -e "\n***** Process Merged annotations for POP: ${pop}, CADD file: $merged_annotation_file"
        python "${launchDir}/bin/process_mergedannotation-final.py" \
             --input_file "${merged_annotation_file}" \
             --output_file "${output_file}" \
             --output_dir "${outputdir}" \
             --population ${pop}
        """
}


workflow PROCESS_MERGED_ANNOTATION {
    take:
    pops_ch
    merged_scaled_annotfile_channel

    main:
    if (params.analysis == "GLOBAL") {
        
        //Channel.value("$params.cohort_name").set { pops_ch }
        //pops_ch.merge(rv_file_ch, vep_loftee_file_ch)

        // Combine the final outputs into a single channel
        combined_process_mergeannot_ch = pops_ch
            .merge(merged_scaled_annotfile_channel)
            .map{ 
                pops, merged_annotfile -> 
                tuple(pops, merged_annotfile)
        }

        combined_process_mergeannot_ch.view({ "SUBWORKFLOW PROCESSED MERGED FINAL WATERSHED ANNOT files combined output: $it" })

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

