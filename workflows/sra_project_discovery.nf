#!/usr/bin/env nextflow

include { RESOLVE_SRA_PROJECT; VALIDATE_SRA_PROJECT } from '../modules/local/sra_project/main'

workflow SRA_PROJECT_DISCOVERY {
    main:
    if (!params.sraProject) {
        error 'SRA_DISCOVERY requires --sra-project'
    }
    if (params.input) {
        error '--input and --sra-project are mutually exclusive'
    }

    ch_resolver = channel.value(
        file("${projectDir}/bin/resolve_sra_project.py", checkIfExists: true)
    )

    RESOLVE_SRA_PROJECT(
        channel.value(params.sraProject),
        channel.value(params.sraPlatforms),
        channel.value(params.sraEmail ?: ''),
        ch_resolver
    )

    VALIDATE_SRA_PROJECT(
        RESOLVE_SRA_PROJECT.out.run_manifest,
        RESOLVE_SRA_PROJECT.out.sample_manifest,
        RESOLVE_SRA_PROJECT.out.exclusions,
        RESOLVE_SRA_PROJECT.out.summary,
        RESOLVE_SRA_PROJECT.out.runinfo,
        ch_resolver
    )

    emit:
    run_manifest    = RESOLVE_SRA_PROJECT.out.run_manifest
    sample_manifest = RESOLVE_SRA_PROJECT.out.sample_manifest
    exclusions      = RESOLVE_SRA_PROJECT.out.exclusions
    summary         = RESOLVE_SRA_PROJECT.out.summary
    runinfo         = RESOLVE_SRA_PROJECT.out.runinfo
    validated       = VALIDATE_SRA_PROJECT.out.sentinel
}
