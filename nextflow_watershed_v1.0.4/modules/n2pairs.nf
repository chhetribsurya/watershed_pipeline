#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process generate_n2_pairs {
    tag { "JOB ${pop}" }
    label 'watershed_py_env'
    publishDir "$params.out_dir/outputs/4_n2pair_files/", mode : 'copy'

    input:
        tuple val(pop), val(rv_file), val(exp_outlier_file)

    output:
        path "n2pair/${pop}/gene-${pop}-rv.N2pairs.tsv", emit: n2pairs
        path "n2pair/${pop}/gene-${pop}-rv*.txt", emit: all

    script:
        def outputdir = "n2pair/${pop}"
        def outfile = "${outputdir}/gene-${pop}-rv.N2pairs.tsv"

        """
        mkdir -p $outputdir
        python "${launchDir}/bin/step15_generateN2pairs-final.py" \
            --rare_variant_file "$rv_file" \
            --exp_outlier_file "$exp_outlier_file" \
            --outfile "$outfile"
        """
}

