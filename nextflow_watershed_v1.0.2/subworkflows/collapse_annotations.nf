#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*include { process_data; call_phylopscore_annotation } from "../modules/call_annotation_process.nf" params(params)*/

process collapse_vep_annotation {
    debug true
    tag { "JOB ${pop}" }
    label 'watershed_r_env'
    publishDir "$params.out_dir/outputs/3_collapse_annotation/", mode : 'copy'

    input:
        tuple val(pop), val(rv_file), val(vep_loftee_file)
        /*tuple val(pop)
        val(rv_file)
        val(vep_loftee_file)*/

    output:
        path "vep/${pop}/gene-${pop}-rv.vep.loftee.collapse.tsv", emit: collapse_out
        path "vep/${pop}/gene-${pop}-rv.vep.loftee.uncollapsed.rvpair.tsv", emit: uncollapse_out
        path "vep/${pop}/annotation_inputfinal/genelevel_output"
        path "vep/${pop}/annotation_inputfinal/variantlevel_output"

    script:
        //def outputdir = "vep/${pop}"
        def genelevel_dir="vep/${pop}/annotation_inputfinal/genelevel_output"
        def variantlevel_dir="vep/${pop}/annotation_inputfinal/variantlevel_output"
        def outfile1 = "vep/${pop}/gene-${pop}-rv.vep.loftee.collapse.tsv"
        def outfile2 = "vep/${pop}/gene-${pop}-rv.vep.loftee.uncollapsed.rvpair.tsv"

        """
        #!/usr/bin/bash
        echo -e "\n***** VEP Annotation collapsing of POP: ${pop}, RV file: $rv_file, VEP Annot file: $vep_loftee_file"
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
    debug true
    tag { "JOB ${pop}" }
    //label 'watershed_cadd_renv'
    publishDir "$params.out_dir/outputs/3_collapse_annotation/", mode : 'copy'

    input:
        tuple val(pop), val(rv_file), val(cadd_base_dir)

    output:
        path "cadd/${pop}/gene-${pop}-rv.CADD.collapse.tsv", emit: collapse_out
        path "cadd/${pop}/gene-${pop}-rv.CADD.uncollapsed.rvpair.tsv", emit: uncollapse_out
        path "cadd/${pop}/annotation_inputfinal/genelevel_output"
        path "cadd/${pop}/annotation_inputfinal/variantlevel_output"

    script:
        //def outputdir = "cadd/${pop}"
        def genelevel_dir="cadd/${pop}/annotation_inputfinal/genelevel_output"
        def variantlevel_dir="cadd/${pop}/annotation_inputfinal/variantlevel_output"
        def outfile1 = "cadd/${pop}/gene-${pop}-rv.CADD.collapse.tsv"
        def outfile2 = "cadd/${pop}/gene-${pop}-rv.CADD.uncollapsed.rvpair.tsv"

        """
        #!/usr/bin/bash
        echo -e "\n***** CADD Annotation collapsing of POP: ${pop}, RV file: $rv_file, CADD file dir: $cadd_base_dir"
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
    debug true
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
        //def outputdir = "ucsc/${pop}"
        def genelevel_dir="ucsc/${pop}/annotation_inputfinal/genelevel_output"
        def variantlevel_dir="ucsc/${pop}/annotation_inputfinal/variantlevel_output"
        def outfile1 = "ucsc/${pop}/gene-${pop}-rv.phyloP100way.collapse.tsv"
        def outfile2 = "ucsc/${pop}/gene-${pop}-rv.phyloP100way.uncollapsed.rvpair.tsv"

        """
        #!/usr/bin/bash
        echo -e "\n***** UCSC Annotation collapsing of POP: ${pop}, RV file: $rv_file, PhyloP100way Annot file: $phyloP100way_file"
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
    debug true
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
        //def outputdir = "gencode/${pop}"
        def genelevel_dir="gencode/${pop}/annotation_inputfinal/genelevel_output"
        def variantlevel_dir="gencode/${pop}/annotation_inputfinal/variantlevel_output"
        def outfile1 = "gencode/${pop}/gene-${pop}-rv.gencode.collapse.tsv"
        def outfile2 = "gencode/${pop}/gene-${pop}-rv.gencode.uncollapsed.rvpair.tsv"

        """
        #!/usr/bin/bash
        echo -e "\n***** GENCODE Annotation collapsing of POP: ${pop}, GENCODE Annot file: $gencode_file"
        Rscript "${launchDir}/bin/step12.4_collapse_gencodeDistFromGene-final.R" \
            "$gencode_file" \
            "$outfile1" \
            "$outfile2" \
            "$genelevel_dir" \
            "$variantlevel_dir"
        """
}

process collapse_afreq_annotation {
    debug true
    tag { "JOB ${pop}" }
    //label 'watershed_r_env'
    publishDir "$params.out_dir/outputs/3_collapse_annotation/", mode : 'copy'

    input:
        tuple val(pop), val(rv_file)

    output:
        path "afreq/${pop}/gene-${pop}-rv.afreq.collapse.tsv", emit: collapse_out
        path "afreq/${pop}/gene-${pop}-rv.afreq.uncollapsed.rvpair.tsv", emit: uncollapse_out
        path "afreq/${pop}/annotation_inputfinal/genelevel_output"
        path "afreq/${pop}/annotation_inputfinal/variantlevel_output"

    script:
        //def outputdir = "afreq/${pop}"
        def genelevel_dir="afreq/${pop}/annotation_inputfinal/genelevel_output"
        def variantlevel_dir="afreq/${pop}/annotation_inputfinal/variantlevel_output"
        def outfile1 = "afreq/${pop}/gene-${pop}-rv.afreq.collapse.tsv"
        def outfile2 = "afreq/${pop}/gene-${pop}-rv.afreq.uncollapsed.rvpair.tsv"

        """
        #!/usr/bin/bash
        echo -e "\n***** Allelic frequency and rarevar number annotation collapsing of POP: ${pop}, Rare variant file: $rv_file"
        Rscript "${launchDir}/bin/step12.5_collapse_aFreq-final.R" \
            "$rv_file" \
            "$outfile1" \
            "$outfile2" \
            "$genelevel_dir" \
            "$variantlevel_dir"
        """
}


workflow COLLAPSE_ANNOTATION {
    take:
    pops_ch
    rv_file_ch
    vep_loftee_file_ch
    cadd_files_dir_ch
    ucsc_phylopfile_ch
    gencode_file_ch 

    main:
    if (params.analysis == "GLOBAL") {
        
        //Channel.value("$params.cohort_name").set { pops_ch }
        //pops_ch.merge(rv_file_ch, vep_loftee_file_ch)
        //Channel.value("$rv_file").set { rv_file_ch }
        //Channel.value("$vep_loftee_file").set {vep_file_ch }

        // Combine the final outputs into a single channel and form tuples
        combined_vep_ch = pops_ch
            .merge(rv_file_ch, vep_loftee_file_ch)
            .map{ pops, rv_file, vep_file -> tuple(pops, rv_file, vep_file)
        }

        combined_vep_ch.view({ "SUBWORKFLOW (COLLAPSE VEP ANNOTATION) combined output: $it" })
        
        // Collapse VEP Annotation
        //collapse_vep_annotation(pops_ch, rv_file_ch, vep_loftee_file_ch)
        collapse_vep_annotation(combined_vep_ch)

        // Combine the final outputs into a single channel and form tuples
        combined_cadd_ch = pops_ch
            .merge(rv_file_ch, cadd_files_dir_ch)
            .map{ pops, rv_file, cadd_dir -> tuple(pops, rv_file, cadd_dir)
        }

        combined_cadd_ch.view({ "SUBWORKFLOW (COLLAPSE CADD ANNOTATION) combined output: $it" })

        // Collapse CADD Annotation
        collapse_cadd_annotation(combined_cadd_ch)


        // Combine the final outputs into a single channel and form tuples
        combined_ucsc_ch = pops_ch
            .merge(rv_file_ch, ucsc_phylopfile_ch)
            .map{ pops, rv_file, phylop_file -> tuple(pops, rv_file, phylop_file)
        }

        combined_ucsc_ch.view({ "SUBWORKFLOW (COLLAPSE UCSC PhyloP ANNOTATION) combined output: $it" })

        // Collapse UCSC phylop score Annotation
        collapse_ucsc_annotation(combined_ucsc_ch)


        // Combine the final outputs into a single channel and form tuples
        combined_gencode_ch = pops_ch
            .merge(gencode_file_ch)
            .map{ pops, gencode_file -> tuple(pops, gencode_file)
        }

        combined_gencode_ch.view({ "SUBWORKFLOW (COLLAPSE GENCODE ANNOTATION) combined output: $it" })

        // Collapse GENCODE annotation
        collapse_gencode_annotation(combined_gencode_ch)


        // Combine the final outputs into a single channel and form tuples
        combined_afreq_ch = pops_ch
            .merge(rv_file_ch)
            .map{ pops, rv_file -> tuple(pops, rv_file)
        }

        combined_afreq_ch.view({ "SUBWORKFLOW (COLLAPSE AFREQ ANNOTATION) combined output: $it" })

        // Collapse Allelic Frequency(afreq) and rare variant number annotation
        collapse_afreq_annotation(combined_afreq_ch)



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
    vep_collapsed_file = collapse_vep_annotation.out.collapse_out
    vep_uncollapsed_file = collapse_vep_annotation.out.uncollapse_out
    cadd_collapsed_file = collapse_cadd_annotation.out.collapse_out
    cadd_uncollapsed_file = collapse_cadd_annotation.out.uncollapse_out
    ucsc_collapsed_file = collapse_ucsc_annotation.out.collapse_out
    ucsc_uncollapsed_file = collapse_ucsc_annotation.out.uncollapse_out
    gencode_collapsed_file = collapse_gencode_annotation.out.collapse_out
    gencode_uncollapsed_file = collapse_gencode_annotation.out.uncollapse_out
    afreq_collapsed_file = collapse_afreq_annotation.out.collapse_out
    afreq_uncollapsed_file = collapse_afreq_annotation.out.uncollapse_out
}


