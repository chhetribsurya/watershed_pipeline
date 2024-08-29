#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
=================================================================================================
    Rare-Variant-Watershed Nextflow Pipeline
=================================================================================================
    @author  : Surya B. Chhetri
    Github   : chhetribsurya@github.com
    Contact  : chhetribsurya@gmail.com
-------------------------------------------------------------------------------------------------
*/


// DEFAULT PRINT WITH PIPELINE RUN
version                 = "1.0.3"
println """\

         R A R E - V A R I A N T - W A T E R S H E D - N F   P I P E L I N E  ~  Version ${version}
         ==========================================================================================
         TISSUE             : ${params.tissue}
         VCF FILE           : ${params.rv_file}
         TPM FILE           : ${params.tpm_infile}
         READCOUNT FILE     : ${params.readcount_infile}
         COVARIATE FILE     : ${params.covariate_infile}
         SUBJID FILE        : ${params.subjids_file}
         OUTPUT DIR         : ${params.out_dir}
         """


log.info '\n'
log.info '==================================================================================================!'
log.info '==================================================================================================!'
log.info '\n'

// Load subworkflow functions
include { CALL_OUTLIERS               } from './subworkflows/call_rv_outliers' params(params)
include { CALL_ANNOTATION             } from './subworkflows/call_annotation' params(params)
include { CALL_RAREVARIANTS           } from './subworkflows/call_rvs' params(params)
include { COLLAPSE_ANNOTATION         } from './subworkflows/collapse_annotations' params(params)
include { CALL_N2_PAIRS               } from './subworkflows/n2pairs' params(params)
include { MERGE_ANNOTATIONS           } from './subworkflows/merge_annotations' params(params)
include { PROCESS_MERGED_ANNOTATION   } from './subworkflows/process_merge_annotations' params(params)
include { SAMPLE_SELECT               } from './subworkflows/samples_to_analyze.nf' params(params)
include { WATERSHED_MODELRUN          } from './subworkflows/watershed_modelrun' params(params)


workflow {
    main:
        //SAMPLE SELECTION FOR ANALYSIS
        Channel.value("$params.cohort_name").set { pops_ch }
        Channel.value("$params.subjids_file").set { subjids_file_ch }

        SAMPLE_SELECT(pops_ch, subjids_file_ch) 
        subjids_channel = SAMPLE_SELECT.out[0]

        //DEFINE RV CHANNELS
        Channel.value("$params.cohort_name").set { pops_ch }
        Channel.fromPath("$params.rv_file").set { rv_file_ch }
        Channel.fromPath("$params.gencode_region_file").set { gencode_region_file_ch }
        Channel.fromPath("$params.gnomad_file").set { gnomad_file_ch }
        Channel.fromPath("$params.gnomad_generegion_filelist").set { gnomad_generegion_filelist_ch }
        Channel.value("$params.ancestry").set { ancestry_channel }

        CALL_RAREVARIANTS(pops_ch, rv_file_ch, gencode_region_file_ch, gnomad_file_ch, \
                        subjids_channel, gnomad_generegion_filelist_ch, ancestry_channel)
        
        //DEFINE RV PARAMETERS
        gene_rarevar_ch  = CALL_RAREVARIANTS.out[10]
        gene_rarevar_ch2 = CALL_RAREVARIANTS.out[9]

        //DEFINE OUTLIER CHANNELS
        CALL_OUTLIERS(pops_ch, subjids_channel, gene_rarevar_ch2)
        
        //DEFINE ANNOTATION CHANNELS
        CALL_ANNOTATION(pops_ch, gene_rarevar_ch, gene_rarevar_ch2)

        // DEFINE ANNOTATION COLLAPSE PARAMETERS
        exp_outlier_file_channel = CALL_OUTLIERS.out[0]
        vep_loftee_file_ch = CALL_ANNOTATION.out[7]
        cadd_files_dir_ch = CALL_ANNOTATION.out[8]
        ucsc_phylop_file_channel = CALL_ANNOTATION.out[0]
        gencode_file_channel = CALL_ANNOTATION.out[1]

        // COLLECT ALL OUTPUTS FROM cadd_files_dir_ch WAITING ON ALL PROCESS
        cadd_files_dir_ch.collect().map { dirs -> dirs.unique() }.set { unique_cadd_files_dir_ch }
        cadd_files_dir_ch.collect().map { dirs -> dirs }.set { nonunique_cadd_files_dir_ch }

        COLLAPSE_ANNOTATION(pops_ch, gene_rarevar_ch2, vep_loftee_file_ch, unique_cadd_files_dir_ch, 
                            ucsc_phylop_file_channel, gencode_file_channel)

        CALL_N2_PAIRS(pops_ch, gene_rarevar_ch2, exp_outlier_file_channel)

        // DEFINE MERGE ANNOTATION PARAMETERS
        cadd_annotfile_ch = COLLAPSE_ANNOTATION.out[2]
        gencode_annotfile_ch = COLLAPSE_ANNOTATION.out[6]
        afreq_annotfile_ch = COLLAPSE_ANNOTATION.out[8]
        ucsc_phylop_annotfile_ch = COLLAPSE_ANNOTATION.out[4]
        vep_annotfile_ch = COLLAPSE_ANNOTATION.out[0]
        outlier_refbasefile_ch = CALL_OUTLIERS.out[3]
        n2pair_file_channel = CALL_N2_PAIRS.out[0]

        MERGE_ANNOTATIONS(pops_ch, cadd_annotfile_ch, gencode_annotfile_ch, afreq_annotfile_ch, ucsc_phylop_annotfile_ch, 
                        vep_annotfile_ch, outlier_refbasefile_ch, n2pair_file_channel)

        // DEFINE ANNOTATION PROCESSING PARAMETERS
        merged_unscaled_annotfile_channel = MERGE_ANNOTATIONS.out[0]
        merged_scaled_annotfile_channel = MERGE_ANNOTATIONS.out[1]

        PROCESS_MERGED_ANNOTATION(pops_ch, merged_scaled_annotfile_channel)

        // DEFINE WATERSHED PARAMETERS
        processed_merged_annotfile_ch = PROCESS_MERGED_ANNOTATION.out[0]
        feature_countfile_channel = PROCESS_MERGED_ANNOTATION.out[1]

        WATERSHED_MODELRUN(pops_ch, processed_merged_annotfile_ch)

}


workflow.onComplete {

    println ( workflow.success ? """
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """ : """
        Failed: ${workflow.errorReport}
        exit status : ${workflow.exitStatus}
        """
    )
}


/*
========================================================================================
    THE END
========================================================================================
*/

