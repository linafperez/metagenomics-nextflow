#!/usr/bin/env nextflow

include { METAGENOMICS } from './workflows/metagenomics'
include { SRA_PROJECT_DISCOVERY } from './workflows/sra_project_discovery'
include { SRA_CHECKPOINT_RECONCILIATION } from './workflows/sra_checkpoint_reconciliation'
include { SRA_SAMPLE_PREPROCESSING } from './workflows/sra_sample_preprocessing'
include { SRA_GLOBAL } from './workflows/sra_global'

workflow {
    def stage = params.executionStage ?: 'auto'
    if (params.input && params.sraProject) {
        error '--input and --sra-project are mutually exclusive'
    }
    if (stage == 'auto') {
        if (params.input && !params.sraProject) {
            METAGENOMICS()
        } else if (params.sraProject) {
            error 'BioProject mode uses disk-safe staged execution; run metagenomics_pipeline.sh --sra-project PRJ...'
        } else {
            error 'exactly one production input is required: --input or --sra-project via the launcher'
        }
    } else if (stage == 'local') {
        if (!params.input || params.sraProject) {
            error 'local stage requires --input and forbids --sra-project'
        }
        METAGENOMICS()
    } else if (stage == 'sra-discovery') {
        if (!params.sraProject || params.input) {
            error 'SRA stages require --sra-project and forbid --input'
        }
        SRA_PROJECT_DISCOVERY()
    } else if (stage == 'sra-checkpoints') {
        if (!params.sraProject || params.input) {
            error 'SRA stages require --sra-project and forbid --input'
        }
        SRA_CHECKPOINT_RECONCILIATION()
    } else if (stage == 'sra-preprocess') {
        if (!params.sraProject || params.input) {
            error 'SRA stages require --sra-project and forbid --input'
        }
        SRA_SAMPLE_PREPROCESSING()
    } else if (stage == 'sra-global') {
        if (!params.sraProject || params.input) {
            error 'SRA stages require --sra-project and forbid --input'
        }
        SRA_GLOBAL()
    } else {
        error "unsupported internal execution stage: ${stage}"
    }
}
