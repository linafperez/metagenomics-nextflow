#!/usr/bin/env nextflow

include { SRA_PROJECT_DISCOVERY } from '../sra_project_discovery/main'
include { SRA_CHECKPOINT_RECONCILIATION } from '../sra_checkpoint_reconciliation/main'
include { SRA_SAMPLE_PREPROCESSING } from '../sra_sample_preprocessing/main'
include { SRA_GLOBAL } from '../sra_global/main'

workflow SRA_INPUT {
    main:
    def stage = params.executionStage?.toString() ?: ''
    if (stage == 'sra-discovery') {
        SRA_PROJECT_DISCOVERY()
    } else if (stage == 'sra-checkpoints') {
        SRA_CHECKPOINT_RECONCILIATION()
    } else if (stage == 'sra-preprocess') {
        SRA_SAMPLE_PREPROCESSING()
    } else if (stage == 'sra-global') {
        SRA_GLOBAL()
    } else {
        error "unsupported internal SRA execution stage: ${stage}"
    }
}
