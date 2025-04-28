#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {
  call_phylopscore_annotation;
  call_gencode_annotation_final;
  chromsplit_rv_file;
  call_cadd_annotation;
  call_vep_annotation_final;
  prepare_vep_parse;
  parse_vep_annotation
} from "../modules/call_annotation" params(params)


workflow CALL_ANNOTATION {
    take:
    pops_annot_ch
    gene_rarevar_ch
    gene_rarevar_ch2


    main:
    if (params.analysis == "GLOBAL") {
        
        combined_ch = pops_annot_ch
            .merge(gene_rarevar_ch)
            .map { pops_annot, gene_rarevar -> tuple(pops_annot, gene_rarevar)
        }

        combined_ch2 = pops_annot_ch
            .merge(gene_rarevar_ch2)
            .map { pops_annot, gene_rarevar -> tuple(pops_annot, gene_rarevar)
        }

        // Phylop annotation
        call_phylopscore_annotation(combined_ch)

        // Gencode annotation
        call_gencode_annotation_final(combined_ch2)

        // Split VCF files by chrom
        chromsplit_rv_file(combined_ch2)

        // Split files output
        vcffiles_ch = chromsplit_rv_file.out.vcfannotations

        vcffiles_ch.flatten()
            .take(params.num_chromosomes)
            .map { vcf -> tuple("$params.cohort_name", vcf) }
            .set { vcf_channel }

        /*Run CADD Annotation chromosome wise as a separate process*/
        call_cadd_annotation(vcf_channel)
        
        /*Run VEP Annotation chromosome wise as a separate process*/
        call_vep_annotation_final(vcf_channel)

        //Prepare VEP parse of VCFs
        Channel.value("$params.cohort_name").set { pop_channel}    
        vep_vcf_ch = call_vep_annotation_final.out.annotations
        vep_vcf_ch
            .collect()
            .set { vcf_list_channel }

        prepare_vep_parse(vcf_list_channel, pop_channel)

        //Parse VEP annotations
        vcf_newlist_channel = prepare_vep_parse.out.vep_vcf_annotations
        vcf_newlist_channel
            .collect()
            .map { files -> tuple("$params.cohort_name", files) }
            .set { combined_vcf_channel }

        parse_vep_annotation(combined_vcf_channel)


       /*------------------------------------------------------------------*/
       /* if (!params.skip_cache) {
            println "BUILDING VEP CACHE ..." 
            build_vep_cache()
            build_vep_cache
                .out
                .cache_files
                .set { cache_dir_ch }
        } else {
            println "SKIPPING VEP CACHE ..."
            Channel
                .fromPath(params.cachedir_vep)
                .set { cache_dir_ch }
        }*/

        /*Channel
            .fromPath('path/to/vcfs/*.vcf')  // Adjust this path to your VCF files
            .flatten()
            .set { vcf_channel }*/

        //vcffiles_ch.flatten().set { vcf_channel }

        /*vcf_channel
           .combine(cache_dir_ch)
           .view { combined -> "VEP CACHE DIR and VCF FILE COMBINED TUPLE: ${combined}" }
           .set { vep_input_ch }*/

        // Map to create tuples instead of combine
        /*vcf_channel
            .map { vcf -> tuple(vcf, cache_dir_ch) }
            .view { combined -> "VEP CACHE DIR and VCF FILE COMBINED TUPLE: ${combined}" }
            .set { vep_input_ch }*/

        //call_vep_annotation_final(vep_input_ch)
        /*call_vep_annotation_final(vcffiles_ch.flatten(), cache_dir_ch )*/

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
    phylop_annot_file  = call_phylopscore_annotation.out.annotations
    gencode_annot_file  = call_gencode_annotation_final.out.annotations
    chromsplit_file  = chromsplit_rv_file.out.annotations
    chromsplit_vcf_file  = chromsplit_rv_file.out.vcfannotations
    cadd_annot_file  = call_cadd_annotation.out.annotations
    vep_annot_file  = call_vep_annotation_final.out.annotations
    vep_annot_bgzipped = prepare_vep_parse.out.vep_vcf_annotations
    vep_annot_parsed_file = parse_vep_annotation.out.parsed_annotations
    cadd_annot_dir  = call_cadd_annotation.out.cadd_dir
    //vep_cache_file = build_vep_cache.out.cache_files
}


