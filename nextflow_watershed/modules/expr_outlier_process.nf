#!/usr/bin/env nextflow

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



