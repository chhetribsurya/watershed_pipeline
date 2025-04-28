#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process collapse_vep_annotation {
    tag { "JOB ${pop}" }
    label 'watershed_r_env'
    publishDir "$params.out_dir/outputs/3_collapse_annotation/", mode : 'copy'

    input:
        tuple val(pop), val(rv_file), val(vep_loftee_file)

    output:
        path "vep/${pop}/gene-${pop}-rv.vep.loftee.collapse.tsv", emit: collapse_out
        path "vep/${pop}/gene-${pop}-rv.vep.loftee.uncollapsed.rvpair.tsv", emit: uncollapse_out
        path "vep/${pop}/annotation_inputfinal/genelevel_output"
        path "vep/${pop}/annotation_inputfinal/variantlevel_output"

    script:
        def genelevel_dir="vep/${pop}/annotation_inputfinal/genelevel_output"
        def variantlevel_dir="vep/${pop}/annotation_inputfinal/variantlevel_output"
        def outfile1 = "vep/${pop}/gene-${pop}-rv.vep.loftee.collapse.tsv"
        def outfile2 = "vep/${pop}/gene-${pop}-rv.vep.loftee.uncollapsed.rvpair.tsv"

        """
        Rscript "${launchDir}/bin/step12.1_collapse_vep_loftee-final.R" \
            "$rv_file" \
            "$outfile1" \
            "$outfile2" \
            "$genelevel_dir" \
            "$variantlevel_dir" \
            "$vep_loftee_file"
        """
}


process collapse_cadd_annotation {
    tag { "JOB ${pop}" }
    publishDir "$params.out_dir/outputs/3_collapse_annotation/", mode : 'copy'

    input:
        tuple val(pop), val(rv_file), val(cadd_base_dir)

    output:
        path "cadd/${pop}/gene-${pop}-rv.CADD.collapse.tsv", emit: collapse_out
        path "cadd/${pop}/gene-${pop}-rv.CADD.uncollapsed.rvpair.tsv", emit: uncollapse_out
        path "cadd/${pop}/annotation_inputfinal/genelevel_output"
        path "cadd/${pop}/annotation_inputfinal/variantlevel_output"

    script:
        def genelevel_dir="cadd/${pop}/annotation_inputfinal/genelevel_output"
        def variantlevel_dir="cadd/${pop}/annotation_inputfinal/variantlevel_output"
        def outfile1 = "cadd/${pop}/gene-${pop}-rv.CADD.collapse.tsv"
        def outfile2 = "cadd/${pop}/gene-${pop}-rv.CADD.uncollapsed.rvpair.tsv"

        """
        Rscript "${launchDir}/bin/step12.2_collapse_cadd_chromwise.R" \
            "$rv_file" \
            "$outfile1" \
            "$outfile2" \
            "$genelevel_dir" \
            "$variantlevel_dir" \
            "$cadd_base_dir" \
            "$pop"
        """
}


process collapse_ucsc_annotation {
    tag { "JOB ${pop}" }
    label 'watershed_r_env'
    publishDir "$params.out_dir/outputs/3_collapse_annotation/", mode : 'copy'

    input:
        tuple val(pop), val(rv_file), val(phyloP100way_file)

    output:
        path "ucsc/${pop}/gene-${pop}-rv.phyloP100way.collapse.tsv", emit: collapse_out
        path "ucsc/${pop}/gene-${pop}-rv.phyloP100way.uncollapsed.rvpair.tsv", emit: uncollapse_out
        path "ucsc/${pop}/annotation_inputfinal/genelevel_output"
        path "ucsc/${pop}/annotation_inputfinal/variantlevel_output"

    script:
        def genelevel_dir="ucsc/${pop}/annotation_inputfinal/genelevel_output"
        def variantlevel_dir="ucsc/${pop}/annotation_inputfinal/variantlevel_output"
        def outfile1 = "ucsc/${pop}/gene-${pop}-rv.phyloP100way.collapse.tsv"
        def outfile2 = "ucsc/${pop}/gene-${pop}-rv.phyloP100way.uncollapsed.rvpair.tsv"

        """
        Rscript "${launchDir}/bin/step12.3_collapse_ucsc-final.R" \
            "$rv_file" \
            "$outfile1" \
            "$outfile2" \
            "$genelevel_dir" \
            "$variantlevel_dir" \
            "$phyloP100way_file"
        """
}


process collapse_gencode_annotation {
    tag { "JOB ${pop}" }
    label 'watershed_r_env'
    publishDir "$params.out_dir/outputs/3_collapse_annotation/", mode : 'copy'

    input:
        tuple val(pop), val(gencode_file)

    output:
        path "gencode/${pop}/gene-${pop}-rv.gencode.collapse.tsv", emit: collapse_out
        path "gencode/${pop}/gene-${pop}-rv.gencode.uncollapsed.rvpair.tsv", emit: uncollapse_out
        path "gencode/${pop}/annotation_inputfinal/genelevel_output"
        path "gencode/${pop}/annotation_inputfinal/variantlevel_output"

    script:
        def genelevel_dir="gencode/${pop}/annotation_inputfinal/genelevel_output"
        def variantlevel_dir="gencode/${pop}/annotation_inputfinal/variantlevel_output"
        def outfile1 = "gencode/${pop}/gene-${pop}-rv.gencode.collapse.tsv"
        def outfile2 = "gencode/${pop}/gene-${pop}-rv.gencode.uncollapsed.rvpair.tsv"

        """
        Rscript "${launchDir}/bin/step12.4_collapse_gencodeDistFromGene-final.R" \
            "$gencode_file" \
            "$outfile1" \
            "$outfile2" \
            "$genelevel_dir" \
            "$variantlevel_dir"
        """
}

process collapse_afreq_annotation {
    tag { "JOB ${pop}" }
    publishDir "$params.out_dir/outputs/3_collapse_annotation/", mode : 'copy'

    input:
        tuple val(pop), val(rv_file)

    output:
        path "afreq/${pop}/gene-${pop}-rv.afreq.collapse.tsv", emit: collapse_out
        path "afreq/${pop}/gene-${pop}-rv.afreq.uncollapsed.rvpair.tsv", emit: uncollapse_out
        path "afreq/${pop}/annotation_inputfinal/genelevel_output"
        path "afreq/${pop}/annotation_inputfinal/variantlevel_output"

    script:
        def genelevel_dir="afreq/${pop}/annotation_inputfinal/genelevel_output"
        def variantlevel_dir="afreq/${pop}/annotation_inputfinal/variantlevel_output"
        def outfile1 = "afreq/${pop}/gene-${pop}-rv.afreq.collapse.tsv"
        def outfile2 = "afreq/${pop}/gene-${pop}-rv.afreq.uncollapsed.rvpair.tsv"

        """
        Rscript "${launchDir}/bin/step12.5_collapse_aFreq-final.R" \
            "$rv_file" \
            "$outfile1" \
            "$outfile2" \
            "$genelevel_dir" \
            "$variantlevel_dir"
        """
}

