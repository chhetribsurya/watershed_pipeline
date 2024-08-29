#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process split_population {
    tag { "JOB ${params.cohort}" }
    label 'watershed_python_env'
    publishDir "$params.out_dir/outputs/split_pops", mode : 'copy'

    input:
    path subjids_file

    output:
    path("*_ids.txt"), emit: pop_ids

    script:
    """
    python ${launchDir}/bin/split_population.py --input $subjids_file --output_dir ${params.out_dir}
    """
}


process process_data {
    debug true
    tag { "JOB ${pop}" }
    label 'watershed_python_env'
    publishDir "$params.out_dir/outputs/1_expression/$pop", mode : 'copy'

    input:
        tuple val(subjids_file), val(pop)

    output:
       path("*.txt"), emit: all
       path("*covariates.txt"), emit: covariates
       path("*reads.txt"), emit: reads
       path("*tpm.txt"), emit: tpms
       val(pop), emit: pops

    script:
        """
        echo -e "\n\n"
        echo "---------------------------------------------------"
        python ${launchDir}/bin/step1_process_data.py --out_dir $params.out_dir --subjids_file $subjids_file --covariate_infile $params.covariate_infile --pseudocount_infile $params.readcount_infile --tpm_infile $params.tpm_infile --population $pop --tissue $params.tissue
        echo ""
        echo "---------------------------------------------------"
        echo -e "\n\n"
        """
}


process normalize_expression {
    tag { "JOB ${pop}" }
    publishDir "$params.out_dir/outputs/1_expression/$pop/PEER", mode : 'copy'

    input:
        path(tpm_file) 
        path(read_file) 
        path(cov_file)
        val(pop)

    output:
        path("*log2.ztrans.txt"), emit: lognormalized
        val(pop), emit: pops_norm

    script:
        """
        Rscript ${launchDir}/bin/step1_preprocess_expression.R --pop $pop --tissue $params.tissue --COV ${cov_file} --TPM_FILE ${tpm_file} --READ_FILE ${read_file} --min_reads 6 --min_tpm 0.1
        """
}


process normalize_expression_global {
    tag { "JOB ${pop}" }
    label 'watershed_r_env'
    publishDir "$params.out_dir/outputs/1_expression/$pop/PEER", mode : 'copy'

    input:
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
        Rscript ${launchDir}/bin/step1_preprocess_expression.R --pop $pop --tissue $params.tissue --COV ${cov_file} --TPM_FILE ${tpm_file} --READ_FILE ${read_file} --min_reads 6 --min_tpm 0.1
        """
}


process compute_peerfactors {
    tag { "PEER_${pop}" }
    label 'peer_docker_env' // You might need to define this label in your config with appropriate resources
    publishDir "$params.out_dir/outputs/1_expression/$pop/PEER", mode : 'copy'

    input:
        path normalized_files
        val pop

    output:
        path("*${pop}_Factors*/factors.tsv"), emit: peer_factors // Adjust as needed based on actual output structure
        val(pop), emit: pops_peer

    script:
        """
        bash ${launchDir}/bin/step2.0_calculate_PEER_factors.sh -p ${normalized_files} -s ${launchDir}/bin
        """
}


process compute_expr_residuals {
    tag { "JOB ${pop}" }
    label 'watershed_r_env'
    publishDir "$params.out_dir/outputs/1_expression/$pop/PEER/residuals", mode : 'copy'

    input:
        tuple val(pop), val(lognorm_expr_file), val(peerfactor_file), val(cov_file)

    output:
        path("${params.tissue}.${pop}_residuals.txt"), emit: residuals
        val(pop), emit: pops_residuals

    script:
        """
        Rscript ${launchDir}/bin/step2.1_calculate_PEER_residuals-final.R --expr "$lognorm_expr_file" --cov "$cov_file" --peer "$peerfactor_file" --out "${params.tissue}.${pop}_residuals.txt" --pcs "${params.genotype_pcs}"
        """
}


process call_expr_outliers {
    tag { "JOB ${pop}" }
    label 'watershed_r_env'
    publishDir "$params.out_dir/outputs/1_expression", mode : 'copy'

    input:
        tuple val(pop), val(residual_expr_file), val(gene_rarevariant_file)

    output:
        path("rv_expoutlier_refbase/${pop}/*.tsv"), emit: outliers
        val(pop), emit: pops_outliers
        path("rv_expoutlier_refbase/${pop}/*expoutlier.Pvalthresbased.reference.tsv"), emit: outlier_pvalthresh
        path("rv_expoutlier_refbase/${pop}/*Pvalthresbased.globalOutliersRemoved*.txt"), emit: outlier_removed_pval
        path("rv_expoutlier_refbase/${pop}/*Zscorethresbased.globalOutliersRemoved*.txt"), emit: outlier_removed_zscore

    script:
        def outputdir = "rv_expoutlier_refbase/${pop}"

        """
        Rscript ${launchDir}/bin/step2_outlier_calling_script-final.R $residual_expr_file $gene_rarevariant_file $outputdir $pop
        """
}


process call_expr_outliers_only {
    tag { "JOB ${pop}" }
    label 'watershed_r_env'
    publishDir "$params.out_dir/outputs/1_expression", mode : 'copy'

    input:
        tuple val(pop), val(residual_expr_file)

    output:
        path("expression_outlier_only/${pop}/*thresbased.globalOutliersRemoved*.txt"), emit: outliers
        val(pop), emit: pops_outliers

    script:
        def outputdir = "expression_outlier_only/${pop}"

        """
        Rscript ${launchDir}/bin/step2_outlier_calling_script-only.R $residual_expr_file $outputdir $pop
        """
}


