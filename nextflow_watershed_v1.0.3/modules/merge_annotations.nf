#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process merge_annotations_and_sortN2pair {
    tag { "JOB ${pop}" }
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
        mkdir -p $annotation_dir
        mkdir -p $outputdir
        cp $cadd_annotfile $gencode_annotfile $afreq_annotfile $ucsc_phylop_annotfile $vep_annotfile $outlier_refbasefile $n2pair_file $outputdir 
        python "${launchDir}/bin/step16_merge_annotations_and_sortN2pair.py" \
            --pop "$pop" \
            --annotation_dir "$annotation_dir" \
            --refbase_file "$outlier_refbasefile" \
            --n2pair_file "$n2pair_file"
        """
}

