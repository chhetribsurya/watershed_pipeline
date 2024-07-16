#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*include { process_data; call_phylopscore_annotation } from "../modules/call_annotation_process.nf" params(params)*/

process merge_annotations_and_sortN2pair {
    debug true
    tag { "JOB ${pop}" }
    //label 'watershed_py_env'
    publishDir "$params.out_dir/outputs/5_merged_annotation/", mode : 'copy'

    input:
        tuple val(pop), val(cadd_annotfile), val(gencode_annotfile), val(afreq_annotfile), val(ucsc_phylop_annotfile), val(vep_annotfile), val(outlier_refbasefile), val(n2pair_file)

    output:
        path "annotation_inputfinal/genelevel_output/${pop}/final-${pop}-rv.mergedannotation.Pvalthresbased.tsv", emit: merged_annotation_only
        path "annotation_inputfinal/genelevel_output/${pop}/final-${pop}-rv.mergedannotation.plusN2pair.Pvalthresbased.tsv", emit: merged_annotation_plusN2pair
        path "annotation_inputfinal/genelevel_output/${pop}/final-${pop}-rv.mergedannotation.plusN2pair.Pvalthresbased.stdscaled.tsv", emit: merged_stdscaled_annotation_plusN2pair
        path "annotation_inputfinal/genelevel_output/${pop}/final-${pop}-rv.mergedannotation.*.tsv", emit: all

    script:
        def annotation_dir = "annotation_inputfinal/genelevel_output"
        def outputdir = "$annotation_dir/${pop}"

        """
        #!/usr/bin/bash
        mkdir -p $annotation_dir
        mkdir -p $outputdir
        echo "Copied annot files to temp annotation dir: $annotation_dir"
        cp $cadd_annotfile $gencode_annotfile $afreq_annotfile $ucsc_phylop_annotfile $vep_annotfile $outlier_refbasefile $n2pair_file $outputdir 
        echo -e "\n***** Merge annotations for POP: ${pop}, CADD file: $cadd_annotfile, GENCODE file: $gencode_annotfile, AFREQ: $afreq_annotfile, PHYLOP file: $ucsc_phylop_annotfile, VEP file: $vep_annotfile, REFBASE file: $outlier_refbasefile, N2PAIR file: $n2pair_file"
        python "${launchDir}/bin/step16_merge_annotations_and_sortN2pair.py" \
            --pop "$pop" \
            --annotation_dir "$annotation_dir" \
            --refbase_file "$outlier_refbasefile" \
            --n2pair_file "$n2pair_file"
        """
}


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
        
        //Channel.value("$params.cohort_name").set { pops_ch }
        //pops_ch.merge(rv_file_ch, vep_loftee_file_ch)

        // Combine the final outputs into a single channel
        combined_mergeannot_ch = pops_ch
            .merge(cadd_annotfile_ch, gencode_annotfile_ch, afreq_annotfile_ch, ucsc_phylop_annotfile_ch, vep_annotfile_ch, outlier_refbasefile_ch, n2pair_file_channel)
            .map{ 
                pops, cadd_annotfile, gencode_annotfile, afreq_annotfile, ucsc_phylop_annotfile, vep_annotfile, outlier_refbasefile, n2pair_file -> 
                tuple(pops, cadd_annotfile, gencode_annotfile, afreq_annotfile, ucsc_phylop_annotfile, vep_annotfile, outlier_refbasefile, n2pair_file)
        }

        combined_mergeannot_ch.view({ "SUBWORKFLOW MERGE VERSION ANNOT files combined output: $it" })

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

