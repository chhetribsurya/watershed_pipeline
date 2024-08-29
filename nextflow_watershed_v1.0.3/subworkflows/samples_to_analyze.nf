#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { identify_samples_to_analyze } from "../modules/samples_to_analyze" params(params)

workflow SAMPLE_SELECT {
    take:
    pops_ch
    subjids_ch    

    main: 
    // Process all individuals together
    input_ch = subjids_ch
        .merge(pops_ch)
        .map{ subjids, pops -> tuple(subjids, pops)
    }

    identify_samples_to_analyze(input_ch)

    /*
    // Collect the outputs from process_data
    tpms_ch = process_data.out.tpms
    reads_ch = process_data.out.reads
    covariates_ch = process_data.out.covariates
    pops_ch = process_data.out.pops

    // Normalize expression process
    normalize_expression(tpms_ch, reads_ch, covariates_ch, pops_ch)
    */

    emit:
    samples_to_analyze  = identify_samples_to_analyze.out.samples
}


