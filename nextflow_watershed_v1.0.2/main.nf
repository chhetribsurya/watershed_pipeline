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

// Show help message
params.help = false //prevents undefined param-warn
if (params.help) {

    log.info '\n'
    log.info '==========================================!'
    log.info 'Watershed\'s NextFlow Pipeline'
    log.info '==========================================!'
    log.info '\n'

    helpMessage()
    
    exit 1
    //exit 0
}

def helpMessage() {
  log.info """
        Usage:
        The typical command for running the pipeline is as follows:
        nextflow run main.nf --query QUERY.vcf --dbDir "processDatabaseDirectory" --dbName "processPrefixName"

        Mandatory arguments:
         --query                        Query file of sequences you wish to process
         --dbDir                        Process database directory (full path required)
         --dbName                       Prefix name of the process database

       Optional arguments:
        --outdir                       Output directory to place final process output
        --outfmt                       Output format ['6']
        --options                      Additional options for process command [-evalue 1e-3]
        --outFileName                  Prefix name for process output [input.processout]
        --threads                      Number of CPUs to use during process job [16]
        --chunkSize                    Number of fasta records to use when splitting the query fasta file
        --app                          Process program to use [processn;processp,tprocessn,processx]
        --help                         This usage statement.
        """
}


// This prevents a warning of undefined parameter
//params.help            = false
params.genome            = false


// DEFAULT PRINT WITH PIPELINE RUN
version                 = "1.0"
println """\

         R A R E - V A R I A N T - W A T E R S H E D - N F   P I P E L I N E  ~  Version ${version}
         ==========================================================================================
         GENOME       : ${params.genome}
         TISSUE       : ${params.tissue}
         VCF FILE     : ${params.inputvcf}
         TPM FILE     : ${params.tpm_infile}
         COHORT       : ${params.variable}
         OUTPUT DIR   : ${params.out_dir}

         """


log.info '\n'
log.info '==================================================================================================!'
log.info '==================================================================================================!'
log.info '\n'


process split_population {
    //cache 'true'
    tag { "JOB ${params.cohort}" }
    label 'watershed_python_env'
    publishDir "$params.out_dir/outputs/split_pops", mode : 'copy'

    input:
    path subjids_file

    output:
    path("*_ids.txt"), emit: pop_ids

    script:
    """
    #python ${params.bin_dir}/split_population.py --input $subjids_file --output_dir ${params.out_dir}
    python ${launchDir}/bin/split_population.py --input $subjids_file --output_dir ${params.out_dir}
    """
}


process process_data {
    //cache 'true'
    debug true
    /*conda "${params.watershed_pyenv}"*/

    //tag { "JOB ${params.variable}" }
    tag { "JOB ${pop}" }
    label 'watershed_python_env'
    publishDir "$params.out_dir/outputs/1_expression/$pop", mode : 'copy'

    input:
        tuple val(subjids_file), val(pop)

    output:
       /*path("*covariates.txt"), optional: true, emit: covariates
       path("*reads.txt"), optional: true, emit: reads
       path("*tpm.txt"), optional: true, emit: tpms*/

       path("*.txt"), emit: all
       path("*covariates.txt"), emit: covariates
       path("*reads.txt"), emit: reads
       path("*tpm.txt"), emit: tpms
       val(pop), emit: pops

    script:
        """
        echo subjidfile: $subjids_file , population : $pop
        #python ${params.bin_dir}/step1_process_data.py --out_dir $params.out_dir --subjids_file $subjids_file --covariate_infile $params.covariate_infile --pseudocount_infile $params.pseudocount_infile --tpm_infile $params.tpm_infile --population $pop --tissue $params.tissue
        python ${launchDir}/bin/step1_process_data.py --out_dir $params.out_dir --subjids_file $subjids_file --covariate_infile $params.covariate_infile --pseudocount_infile $params.readcount_infile --tpm_infile $params.tpm_infile --population $pop --tissue $params.tissue
        """
}


process normalize_expression {
    //cache 'true'
    debug true
    tag { "JOB ${pop}" }
    label 'watershed_r_env'
    publishDir "$params.out_dir/outputs/1_expression/$pop/PEER", mode : 'copy'

    input:
        //tuple path(tpm), path(readcount), path(covariate), val(pop)
        path(tpm_file) 
        path(read_file) 
        path(cov_file)
        val(pop)

    output:
        path("*log2.ztrans.txt"), emit: lognormalized
        val(pop), emit: pops_norm

    script:
        """
        echo "tpm file: $tpm_file, read file: $read_file, cov file: $cov_file, tissue: $params.tissue, population : $pop"
        Rscript ${launchDir}/bin/step1_preprocess_expression.R --pop $pop --tissue $params.tissue --COV ${cov_file} --TPM_FILE ${tpm_file} --READ_FILE ${read_file} --min_reads 6 --min_tpm 0.1
        """
}


process normalize_expression_global {
    //cache 'true'
    debug true
    tag { "JOB ${pop}" }
    label 'watershed_r_env'
    publishDir "$params.out_dir/outputs/1_expression/$pop/PEER", mode : 'copy'

    input:
        //tuple path(tpm), path(readcount), path(covariate), val(pop)
        path(tpm_file) 
        path(read_file) 
        path(cov_file)
        val(pop)
        val(waitstring)

    output:
        path("*log2.ztrans.txt"), emit: lognormalized
        val(pop), emit: pops_norm

    script:
        """
        echo "tpm file: $tpm_file, read file: $read_file, cov file: $cov_file, tissue: $params.tissue, population : $pop"
        Rscript ${launchDir}/bin/step1_preprocess_expression.R --pop $pop --tissue $params.tissue --COV ${cov_file} --TPM_FILE ${tpm_file} --READ_FILE ${read_file} --min_reads 6 --min_tpm 0.1
        """
}


process compute_peerfactors {
    //cache 'true'
    debug true
    tag { "PEER_${pop}" }
    label 'peer_docker_env' // You might need to define this label in your config with appropriate resources
    publishDir "$params.out_dir/outputs/1_expression/$pop/PEER", mode : 'copy'

    input:
        path normalized_files
        val pop

    output:
        path("*${pop}_Factors*/factors.tsv"), emit: peer_factors // Adjust as needed based on actual output structure
        //path("*Factors*/factors.tsv"), emit: peer_factors // Adjust as needed based on actual output structure
        val(pop), emit: pops_peer

    script:
        """
        echo "Compute peerfactors: ${normalized_files} for pop $pop"
        bash ${launchDir}/bin/step2.0_calculate_PEER_factors.sh -p ${normalized_files} -s ${launchDir}/bin
        """
}


process compute_expr_residuals {
    //cache 'true'
    debug true
    tag { "JOB ${pop}" }
    label 'watershed_r_env'
    publishDir "$params.out_dir/outputs/1_expression/$pop/PEER/residuals", mode : 'copy'

    input:
        tuple val(pop), val(lognorm_expr_file), val(peerfactor_file), val(cov_file)
        /*path(lognorm_expr_file) 
        path(cov_file)
        path(peerfactor_file) 
        val(pop)*/

    output:
        path("${params.tissue}.${pop}_residuals.txt"), emit: residuals
        val(pop), emit: pops_residuals

    script:
        """
        echo "lognorm expr file: $lognorm_expr_file, covariate file: $cov_file, peer file: $peerfactor_file, tissue: $params.tissue, population : $pop"
        Rscript ${launchDir}/bin/step2.1_calculate_PEER_residuals-final.R --expr "$lognorm_expr_file" --cov "$cov_file" --peer "$peerfactor_file" --out "${params.tissue}.${pop}_residuals.txt" --pcs "${params.genotype_pcs}"
        """
}


process call_expr_outliers {
    //cache 'true'
    debug true
    tag { "JOB ${pop}" }
    label 'watershed_r_env'
    publishDir "$params.out_dir/outputs/1_expression", mode : 'copy'

    input:
        tuple val(pop), val(residual_expr_file), val(gene_rarevariant_file)

    output:
        path("rv_expoutlier_refbase/${pop}/*.tsv"), emit: outliers
        val(pop), emit: pops_outliers

    script:
        //def outfile = "variant_effect_${vcf_file.baseName}.txt"
        def outputdir = "rv_expoutlier_refbase/${pop}"

        """
        echo "residual expr file: $residual_expr_file, rarevariant file: $gene_rarevariant_file, population : $pop"
        Rscript ${launchDir}/bin/step2_outlier_calling_script-final.R $residual_expr_file $gene_rarevariant_file $outputdir $pop
        """
}


process call_phylopscore_annotation {
    debug true
    tag { "JOB ${pop}" }
    label 'watershed_python_env'
    publishDir "$params.out_dir/outputs/2_annotation", mode : 'copy'

    input:
        tuple val(pop), val(gene_rarevariant_file)

    output:
        path("ucsc/${pop}/*.bed"), emit: annotations
        val(pop), emit: pops_annotations

    script:
        def outputdir = "ucsc/${pop}"

        """
        echo "phylop file: ${params.phylopannot_file}, rarevariant file: $rarevariant_file, population : $pop"
        echo "querying phylop 100 way scores for population: $pop"
        python ${launchDir}/step9.2_phylop100way_query.py --phylop "${params.phylopannot_file}" --rarevariant $rarevariantsbed --pop $pop --outdir $outputdir 
        """
}


/*include { process_data; normalize_expression; compute_peerfactors; compute_expr_residuals;call_expr_outliers } from "./modules/expr_outlier_process" params(params) */
/*include {CALL_OUTLIERS; SIZES} from './subworkflows/call_rv_outliers' params(params)*/


include { CALL_OUTLIERS               } from './subworkflows/call_rv_outliers' params(params)
include { CALL_ANNOTATION             } from './subworkflows/call_annotation' params(params)
include { CALL_RAREVARIANTS           } from './subworkflows/call_rvs' params(params)
include { COLLAPSE_ANNOTATION         } from './subworkflows/collapse_annotations' params(params)
include { CALL_N2_PAIRS               } from './subworkflows/n2pairs' params(params)
include { MERGE_ANNOTATIONS           } from './subworkflows/merge_annotations' params(params)
include { PROCESS_MERGED_ANNOTATION   } from './subworkflows/process_merge_annotations' params(params)
include { WATERSHED_MODELRUN          } from './subworkflows/watershed_modelrun' params(params)


workflow {
    main:
        //DEFINE RV CHANNELS
        Channel.value("$params.cohort_name").set { pops_ch }
        Channel.value("$params.subjids_file").set { subjids_channel }
        Channel.fromPath("$params.rv_file").set { rv_file_ch }
        Channel.fromPath("$params.gencode_region_file").set { gencode_region_file_ch }
        Channel.fromPath("$params.gnomad_file").set { gnomad_file_ch }
        Channel.fromPath("$params.gnomad_generegion_filelist").set { gnomad_generegion_filelist_ch }
        Channel.value("$params.ancestry").set { ancestry_channel }


        CALL_RAREVARIANTS(pops_ch, rv_file_ch, gencode_region_file_ch, gnomad_file_ch, \
                        subjids_channel, gnomad_generegion_filelist_ch, ancestry_channel)
        
        //DEFINE RV PARAMETERS
        gene_rarevar_ch  = CALL_RAREVARIANTS.out[10] //annotation ".rv.bed" file
        gene_rarevar_ch2 = CALL_RAREVARIANTS.out[9] //annotation ".rv.txt" file, equiv to rv_file_ch

        //DEFINE OUTLIER CHANNELS
        Channel.value("$params.cohort_name").set { pops_ch }
        Channel.value("$params.subjids_file").set { subjids_channel }

        CALL_OUTLIERS(pops_ch, subjids_channel, gene_rarevar_ch2)
        
        //DEFINE ANNOTATION CHANNELS
        Channel.value("$params.cohort_name").set { pops_ch }

        CALL_ANNOTATION(pops_ch, gene_rarevar_ch, gene_rarevar_ch2)

        // DEFINE ANNOTATION COLLAPSE PARAMETERS
        exp_outlier_file_channel = CALL_OUTLIERS.out[0]
        Channel.value("$params.cohort_name").set { pops_ch }

        vep_loftee_file_ch = CALL_ANNOTATION.out[7]
        cadd_files_dir_ch = CALL_ANNOTATION.out[8]
        ucsc_phylop_file_channel = CALL_ANNOTATION.out[0]
        gencode_file_channel = CALL_ANNOTATION.out[1]

        // COLLECT ALL OUTPUTS FROM cadd_files_dir_ch WAITING ON ALL PROCESS
        cadd_files_dir_ch.collect().map { dirs -> dirs.unique() }.set { unique_cadd_files_dir_ch }
        cadd_files_dir_ch.collect().map { dirs -> dirs }.set { nonunique_cadd_files_dir_ch }

        COLLAPSE_ANNOTATION(pops_ch, gene_rarevar_ch2, vep_loftee_file_ch, unique_cadd_files_dir_ch, ucsc_phylop_file_channel, gencode_file_channel)

        CALL_N2_PAIRS(pops_ch, gene_rarevar_ch2, exp_outlier_file_channel)

        // DEFINE MERGE ANNOTATION PARAMETERS
        cadd_annotfile_ch = COLLAPSE_ANNOTATION.out[2]
        gencode_annotfile_ch = COLLAPSE_ANNOTATION.out[6]
        afreq_annotfile_ch = COLLAPSE_ANNOTATION.out[8]
        ucsc_phylop_annotfile_ch = COLLAPSE_ANNOTATION.out[4]
        vep_annotfile_ch = COLLAPSE_ANNOTATION.out[0]
        outlier_refbasefile_ch = CALL_OUTLIERS.out[3]
        n2pair_file_channel = CALL_N2_PAIRS.out[0]

        MERGE_ANNOTATIONS(pops_ch, cadd_annotfile_ch, gencode_annotfile_ch, afreq_annotfile_ch, ucsc_phylop_annotfile_ch, vep_annotfile_ch, outlier_refbasefile_ch, n2pair_file_channel)

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

