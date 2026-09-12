#!/usr/bin/env nextflow

include { RESOLVE_SRA_PROJECT; VALIDATE_SRA_PROJECT } from '../../../modules/local/sra_project/main'

workflow SRA_PROJECT_DISCOVERY {
    main:
    if (!params.sraProject || !params.sraSamples) {
        error 'SRA discovery requires --sra-project and --sra-samples'
    }

    def group_column = params.groupColumn?.toString()?.trim() ?: ''
    ch_resolver = channel.value(
        file("${projectDir}/bin/resolve_sra_project.py", checkIfExists: true)
    )
    ch_selection = channel.value(file(params.sraSamples, checkIfExists: true))

    RESOLVE_SRA_PROJECT(
        channel.value(params.sraProject),
        channel.value(params.sraPlatforms),
        channel.value(params.sraEmail ?: ''),
        ch_selection,
        channel.value(group_column),
        ch_resolver
    )

    VALIDATE_SRA_PROJECT(
        RESOLVE_SRA_PROJECT.out.run_manifest,
        RESOLVE_SRA_PROJECT.out.sample_manifest,
        RESOLVE_SRA_PROJECT.out.exclusions,
        RESOLVE_SRA_PROJECT.out.summary,
        RESOLVE_SRA_PROJECT.out.runinfo,
        RESOLVE_SRA_PROJECT.out.requested,
        RESOLVE_SRA_PROJECT.out.metadata,
        ch_selection,
        channel.value(group_column),
        ch_resolver
    )

    emit:
    run_manifest    = RESOLVE_SRA_PROJECT.out.run_manifest
    sample_manifest = RESOLVE_SRA_PROJECT.out.sample_manifest
    requested       = RESOLVE_SRA_PROJECT.out.requested
    sample_metadata = RESOLVE_SRA_PROJECT.out.metadata
    exclusions      = RESOLVE_SRA_PROJECT.out.exclusions
    summary         = RESOLVE_SRA_PROJECT.out.summary
    runinfo         = RESOLVE_SRA_PROJECT.out.runinfo
    validated       = VALIDATE_SRA_PROJECT.out.sentinel
}
