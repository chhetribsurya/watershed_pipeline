#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process process_merged_annotations {
    tag { "JOB ${pop}" }
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
        mkdir -p $outputdir
        python "${launchDir}/bin/process_mergedannotation-final.py" \
             --input_file "${merged_annotation_file}" \
             --output_file "${output_file}" \
             --output_dir "${outputdir}" \
             --population ${pop}
        """
}

