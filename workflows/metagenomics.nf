#!/usr/bin/env nextflow

include { LOCAL_INPUT } from '../subworkflows/local/local_input/main'
include { SRA_INPUT } from '../subworkflows/local/sra_input/main'

workflow METAGENOMICS {
    main:
    def stage = params.executionStage ?: 'auto'
    def has_local = params.input != null && !params.input.toString().trim().isEmpty()
    def has_project = params.sraProject != null && !params.sraProject.toString().trim().isEmpty()
    def has_selection = params.sraSamples != null && !params.sraSamples.toString().trim().isEmpty()
    def group_column = params.groupColumn?.toString()?.trim() ?: ''

    if (group_column && !(group_column ==~ /^[A-Za-z_][A-Za-z0-9_.-]*$/)) {
        error '--group-column contains unsupported characters'
    }

    if (has_local && (has_project || has_selection)) {
        error '--input and --sra-project/--sra-samples are mutually exclusive'
    }
    if (has_project != has_selection) {
        error 'SRA mode requires both --sra-project and --sra-samples'
    }

    if (stage == 'auto') {
        if (has_local) {
            LOCAL_INPUT()
        } else if (has_project) {
            error 'SRA mode uses disk-safe staged execution; run metagenomics_pipeline.sh with --sra-project and --sra-samples'
        } else {
            error 'exactly one production input is required: --input or --sra-project plus --sra-samples'
        }
    } else if (stage == 'local') {
        if (!has_local || has_project || has_selection) {
            error 'local stage requires --input and forbids SRA inputs'
        }
        LOCAL_INPUT()
    } else if (stage in ['sra-discovery', 'sra-checkpoints', 'sra-preprocess', 'sra-global']) {
        if (!has_project || !has_selection || has_local) {
            error 'SRA stages require --sra-project plus --sra-samples and forbid --input'
        }
        SRA_INPUT()
    } else {
        error "unsupported internal execution stage: ${stage}"
    }
}
