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
    publishDir "$params.out_dir/outputs/1_expression/$pop", mode : 'copy'

    input:
        tuple val(pop), val(residual_expr_file), val(gene_rarevariant_file)

    output:
        path("*.txt"), emit: outliers
        val(pop), emit: pops_outliers

    script:
        //def outfile = "variant_effect_${vcf_file.baseName}.txt"
        def outputdir = "rvexpression_${pop}"

        """
        echo "residual expr file: $residual_expr_file, rarevariant file: $gene_rarevariant_file, population : $pop"
        Rscript ${launchDir}/bin/step2_outlier_calling_script-final.R $residual_expr_file $gene_rarevariant_file $outputdir $pop
        """
}


workflow {

    if (params.analysis == "GLOBAL") {
        
        // Process all individuals together
        input_ch = tuple(params.subjids_file, params.cohort_name)
        process_data(input_ch) 

        // specific emit based output
        process_data.out.all.view({"\nHere's the list of ALL processed output file: $it \n"}) //outputs all (everything)
        
        // Collect the outputs from process_data
        tpms_ch = process_data.out.tpms
        reads_ch = process_data.out.reads
        covariates_ch = process_data.out.covariates
        pops_ch = process_data.out.pops

        // Normalize expression process
        normalize_expression(tpms_ch, reads_ch, covariates_ch, pops_ch)

        // Compute PEER factors
        lognorm_expr_ch = normalize_expression.out.lognormalized
        pops_norm_ch = normalize_expression.out.pops_norm

        compute_peerfactors(lognorm_expr_ch, pops_norm_ch)
        compute_peerfactors.out.peer_factors.view({"Here's the PEER factor file: $it"})
        compute_peerfactors.out.pops_peer.view({"Here's the PEER factor pop: $it"})

        // Compute expression residuals
        peer_factors_ch = compute_peerfactors.out.peer_factors
        pops_peer_ch = compute_peerfactors.out.pops_peer


        // Combine the final outputs into a single channel and form tuples
        combined_ch = pops_norm_ch
            .merge(lognorm_expr_ch, peer_factors_ch, covariates_ch) { pops_norm, lognorm_expr, peer_factors, covariates ->
            tuple(pops_norm, lognorm_expr, peer_factors, covariates)
        }
        combined_ch.view({ pops_norm, lognorm_expr, peer_factors, covariates -> "Final combined output: ($pops_norm, $lognorm_expr, $peer_factors, $covariates)" })
        
        compute_expr_residuals(combined_ch)
        compute_expr_residuals.out.residuals.view({"Here's the residual factor file: $it"})
        compute_expr_residuals.out.pops_residuals.view({"Here's the residual factor pop: $it"})

        // Call expression outliers
        residual_expr_ch = compute_expr_residuals.out.residuals
        pops_peer_ch = compute_expr_residuals.out.pops_residuals

        // Combine the final outputs into a single channel and form tuples
        Channel.fromPath($params.gene_rarevar_file).set(gene_rarevar_ch)

        /*Channel.fromPath('./data/*.vcf')
            | map { tuple( it.baseName, it ) }
            | set { sample_vcf_files }*/

        combined_ch = pops_peer_ch
            .merge(residual_expr_ch, gene_rarevar_ch) { pops_norm, residual_expr, gene_rarevar ->
            tuple(pops_norm, residual_expr, gene_rarevar)
        }
        combined_ch.view({ pops_norm, residual_expr, gene_rarevar  -> "Final outlier call step combined output: ($pops_norm, $residual_expr, $gene_rarevar)" })
        
        call_expr_outliers(combined_ch)
        call_expr_outliers.out.outliers.view({"Here's the outlier file list: $it"})
        call_expr_outliers.out.pops_outliers.view({"Here's the outlier pop: $it"})


    } else if (params.analysis == "STRATIFIED") {

        // Split the population file into individual population files and process each one separately
        split_ch = split_population(params.subjids_file).flatten()

        // Create channel
        pop_ch = split_ch.map { file -> tuple(file, file.baseName.split('_')[0]) } // extract population label
        pop_ch.view({"\nHere's the tuple based pop_ch: $it \n"})

        // Process population data
        process_data(pop_ch) 

        /* Collect the outputs from process_data */
        tpms_ch = process_data.out.tpms
        reads_ch = process_data.out.reads
        covariates_ch = process_data.out.covariates
        pops_ch = process_data.out.pops

        // Normalize expression process
        normalize_expression(tpms_ch, reads_ch, covariates_ch, pops_ch)

        lognorm_expr_ch = normalize_expression.out.lognormalized
        pops_norm_ch = normalize_expression.out.pops_norm

        // Compute PEER factors
        compute_peerfactors(lognorm_expr_ch, pops_norm_ch)

        peer_factors_ch = compute_peerfactors.out.peer_factors
        pops_peer_ch = compute_peerfactors.out.pops_peer

        compute_peerfactors.out.peer_factors.view({"Here's the PEER factor file: $it"})
        compute_peerfactors.out.pops_peer.view({"Here's the PEER factor pop: $it"})


        // Keying lognorm_expr_ch using pops_norm_ch
        lognorm_expr_ch = normalize_expression.out.lognormalized
        pops_norm_ch = normalize_expression.out.pops_norm

        lognorm_expr_keyed = lognorm_expr_ch
            .merge(pops_norm_ch) { lognorm_expr, pops_norm -> tuple(pops_norm, lognorm_expr) }
        lognorm_expr_keyed.view({ key, lognorm_expr -> "Keyed lognorm_expr_ch: ($key, $lognorm_expr)" })

        // Keying covariates_ch using pops_ch
        covariates_ch = process_data.out.covariates
        pops_ch = process_data.out.pops

        covariates_keyed = covariates_ch
            .merge(pops_ch) { covariates, pops -> tuple(pops, covariates) }
        covariates_keyed.view({ key, covariates -> "Keyed covariates_ch: ($key, $covariates)" })

        // Keying peer_factors_ch using pops_peer_ch
        peer_factors_ch = compute_peerfactors.out.peer_factors
        pops_peer_ch = compute_peerfactors.out.pops_peer

        peer_factors_keyed = peer_factors_ch
            .merge(pops_peer_ch) { peer_factors, pops_peer -> tuple(pops_peer, peer_factors) }
        peer_factors_keyed.view({ key, peer_factors -> "Keyed peer_factors_ch: ($key, $peer_factors)" })

        // Joining channels based on the key
        combined_ch = lognorm_expr_keyed
            .join(peer_factors_keyed)
            .join(covariates_keyed)
            .map { key, lognorm_expr, peer_factors, covariates ->
                tuple(key, lognorm_expr, peer_factors, covariates)
            }
            .filter { key, lognorm_expr, peer_factors, covariates -> key == "EUR" || key == "AFR" }
            .view({ key, lognorm_expr, peer_factors, covariates -> "Final combined tuple: ($key, $lognorm_expr, $peer_factors, $covariates)" })


        // Compute residuals using the synchronized combined channel
        compute_expr_residuals(combined_ch)
        compute_expr_residuals.out.residuals.view({"Here's the residual factor file: $it"})
        compute_expr_residuals.out.pops_residuals.view({"Here's the residual factor pop: $it"})
        

        // Call expression outliers
        residual_expr_ch = compute_expr_residuals.out.residuals
        pops_expr_ch = compute_expr_residuals.out.pops_residuals
        
        // Keying residual_expr_ch using pops_expr_ch
        residual_expr_keyed = residual_expr_ch
            .merge(pops_expr_ch) { residual_expr, pops_expr -> tuple(pops_expr, residual_expr) }
        residual_expr_keyed.view({ key, residual_expr -> "Keyed residual_expr_ch: ($key, $residual_expr)" })

        // Combine the final outputs into a single channel and form tuples
        // Channel.fromPath($params.gene_rarevar_file).set(gene_rarevar_ch)

        // Keying gene rare variant files with pop-names
        Channel.fromPath('./data/*.vcf')
            | map { tuple( it.baseName.split('-')[1], it ) } // extract population label
            | set { gene_rarevar_keyed }

        // Joining channels based on the key
        combined_ch = residual_expr_keyed
            .join(gene_rarevar_keyed)
            .map { key, residual_expr, gene_rarevar ->
                tuple(key, residual_expr, gene_rarevar)
            }
            .filter { key, residual_expr, gene_rarevar -> key == "EUR" || key == "AFR" }
            .view({ key, residual_expr, gene_rarevar -> "Final combined tuple: ($key, $residual_expr, $gene_rarevar)" })
        
        // Call expression outliers using the synchronized combined channel
        call_expr_outliers(combined_ch)
        call_expr_outliers.out.outliers.view({"Here's the outlier file list: $it"})
        call_expr_outliers.out.pops_outliers.view({"Here's the outlier pop: $it"})

    }

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

