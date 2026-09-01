#!/usr/bin/env nextflow

include { CHECK_SRA_CHECKPOINTS } from '../modules/local/sra_preprocessing/main'

workflow SRA_CHECKPOINT_RECONCILIATION {
    main:
    if (!params.sraManifest || !params.sraCheckpointDir) {
        error 'SRA checkpoint reconciliation requires --sraManifest and --sraCheckpointDir'
    }
    CHECK_SRA_CHECKPOINTS(
        channel.value(file(params.sraManifest, checkIfExists: true)),
        channel.value(file("${projectDir}/bin/manage_sra_checkpoints.py", checkIfExists: true)),
        channel.value(params.sraCheckpointDir),
        channel.value(params.sraRequireComplete.toString().toBoolean())
    )

    emit:
    manifest = CHECK_SRA_CHECKPOINTS.out.manifest
    pending  = CHECK_SRA_CHECKPOINTS.out.pending
    status   = CHECK_SRA_CHECKPOINTS.out.status
}
