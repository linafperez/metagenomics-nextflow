#!/usr/bin/env nextflow

include { LOCAL_INPUT } from '../subworkflows/local/local_input/main'
include { SRA_INPUT } from '../subworkflows/local/sra_input/main'

workflow METAGENOMICS {
    main:
    def has_local = params.input != null && !params.input.toString().trim().isEmpty()
    def has_project = params.sraProject != null && !params.sraProject.toString().trim().isEmpty()
    def has_selection = params.sraSamples != null && !params.sraSamples.toString().trim().isEmpty()
    def group_column = params.groupColumn?.toString()?.trim() ?: ''
    def spades_mode = params.spadesCoassemblyMode?.toString()

    if (!(spades_mode in ['global', 'condition'])) error '--spadesCoassemblyMode must be global or condition'
    if (spades_mode == 'condition' && !group_column) error '--spadesCoassemblyMode condition requires --groupColumn'
    if (group_column && !(group_column ==~ /^[A-Za-z_][A-Za-z0-9_.-]*$/)) error 'Invalid --groupColumn'
    if (has_local && (has_project || has_selection)) error '--input and --sraProject/--sraSamples are mutually exclusive'
    if (has_project != has_selection) error 'SRA mode requires both --sraProject and --sraSamples'
    if (!params.host_bowtie2_index) error '--host_bowtie2_index is required'
    if (!(params.storage_max_jobs.toString() ==~ /^[1-9][0-9]*$/)) error '--storage_max_jobs must be positive'
    if (params.publish_dir_mode == 'move') error '--publish_dir_mode move removes cached outputs; use copy or link to preserve -resume'

    if (has_local) {
        LOCAL_INPUT()
    } else if (has_project) {
        if (!(params.sraProject.toString() ==~ /^PRJ(?:NA|EB|DB)[0-9]+$/)) error 'Invalid --sraProject BioProject accession'
        if (!(params.sraPlatforms.toString() ==~ /^[A-Za-z0-9_.-]+(?:,[A-Za-z0-9_.-]+)*$/)) error 'Invalid --sraPlatforms allowlist'
        if (!(params.sraMaxSize.toString() ==~ /^(?:u|[1-9][0-9]*(?:[KMGT]B?)?)$/)) error 'Invalid --sraMaxSize prefetch limit'
        if (params.sraEmail && !(params.sraEmail.toString() ==~ /^[A-Za-z0-9._%+\-]+@[A-Za-z0-9.\-]+$/)) error 'Invalid --sraEmail'
        if (!(params.sraDownloadMaxForks.toString() ==~ /^[1-9][0-9]*$/)) error '--sraDownloadMaxForks must be positive'
        SRA_INPUT()
    } else {
        error 'Exactly one input is required: --input or --sraProject plus --sraSamples'
    }
}
