#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {
  collapse_vep_annotation;
  collapse_cadd_annotation;
  collapse_ucsc_annotation;
  collapse_gencode_annotation;
  collapse_afreq_annotation
} from "../modules/collapse_annotations" params(params)

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
        
        // Combine the final outputs into a single channel and form tuples
        combined_vep_ch = pops_ch
            .merge(rv_file_ch, vep_loftee_file_ch)
            .map{ pops, rv_file, vep_file -> tuple(pops, rv_file, vep_file)
        }

        // Collapse VEP Annotation
        //collapse_vep_annotation(pops_ch, rv_file_ch, vep_loftee_file_ch)
        collapse_vep_annotation(combined_vep_ch)

        // Combine the final outputs into a single channel and form tuples
        combined_cadd_ch = pops_ch
            .merge(rv_file_ch, cadd_files_dir_ch)
            .map{ pops, rv_file, cadd_dir -> tuple(pops, rv_file, cadd_dir)
        }

        // Collapse CADD Annotation
        collapse_cadd_annotation(combined_cadd_ch)


        // Combine the final outputs into a single channel and form tuples
        combined_ucsc_ch = pops_ch
            .merge(rv_file_ch, ucsc_phylopfile_ch)
            .map{ pops, rv_file, phylop_file -> tuple(pops, rv_file, phylop_file)
        }

        // Collapse UCSC phylop score Annotation
        collapse_ucsc_annotation(combined_ucsc_ch)


        // Combine the final outputs into a single channel and form tuples
        combined_gencode_ch = pops_ch
            .merge(gencode_file_ch)
            .map{ pops, gencode_file -> tuple(pops, gencode_file)
        }

        // Collapse GENCODE annotation
        collapse_gencode_annotation(combined_gencode_ch)


        // Combine the final outputs into a single channel and form tuples
        combined_afreq_ch = pops_ch
            .merge(rv_file_ch)
            .map{ pops, rv_file -> tuple(pops, rv_file)
        }

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


