#!/usr/bin/env nextflow

include { SPADES } from '../../../../modules/core/spades/main'
include { METAQUAST as METAQUAST_SPADES } from '../../../../modules/core/metaquast/main'
include { MERGE_SPADES_ASSEMBLIES } from '../../../../modules/local/merge_spades_assemblies/main'

workflow SPADES_ASSEMBLY {
    take:
    ch_filtered_reads

    main:
    if (params.spadesCoassemblyMode == 'condition') {
        ch_coassembly_input = ch_filtered_reads
            .toList()
            .flatMap { records ->
                SpadesCoassembly.conditionInputs(records).collect { record -> tuple(record[0], record[1]) }
            }
    } else {
        ch_coassembly_input = ch_filtered_reads
            .collect(flat: false)
            .map { records ->
                if (!records) {
                    error "SPAdes coassembly requires at least one paired-end sample"
                }

                def ordered_records = records.toList().sort { left, right ->
                    left[0]['id'].toString() <=> right[0]['id'].toString()
                }

                ordered_records.each { record ->
                    if (record.size() != 2 || !record[0]['id'] || record[1].size() != 2) {
                        error "SPAdes coassembly expects tuples of [meta, [read_1, read_2]]"
                    }
                }

                def sample_ids = ordered_records.collect { record -> record[0]['id'].toString() }
                def reads      = ordered_records.collectMany { record -> record[1].toList() }
                def meta       = [
                    id         : 'spades_coassembly',
                    assembler  : 'spades',
                    branch     : 'spades',
                    sample_ids : sample_ids
                ]

                tuple(meta, reads)
            }
    }

    SPADES(ch_coassembly_input)
    ch_contig_provenance = channel.empty()
    ch_assembly_manifest = channel.empty()
    ch_merge_versions = channel.empty()
    if (params.spadesCoassemblyMode == 'condition') {
        // Collect once after all independent assemblies finish; pair metadata
        // and paths in one ordered record list so arrival order cannot mix them.
        ch_merge_input = SPADES.out.contigs
            .collect(flat: false)
            .map { records ->
                def ordered = records.toList().sort { left, right -> left[0].group <=> right[0].group }
                def meta = [
                    id: 'spades_coassembly', assembler: 'spades', branch: 'spades',
                    coassembly_mode: 'condition', group_column: params.groupColumn,
                    groups: ordered.collect { record -> record[0].group },
                    sample_ids: ordered.collectMany { record -> record[0].sample_ids }.sort()
                ]
                tuple(meta, ordered.collect { record -> record[0] }, ordered.collect { record -> record[1] })
            }
        MERGE_SPADES_ASSEMBLIES(
            ch_merge_input,
            channel.value(file("${projectDir}/bin/spades_coassembly.py", checkIfExists: true))
        )
        ch_assembly = MERGE_SPADES_ASSEMBLIES.out.contigs
        ch_contig_provenance = MERGE_SPADES_ASSEMBLIES.out.provenance
        ch_assembly_manifest = MERGE_SPADES_ASSEMBLIES.out.manifest
        ch_merge_versions = MERGE_SPADES_ASSEMBLIES.out.versions
    } else {
        ch_assembly = SPADES.out.contigs
    }
    METAQUAST_SPADES(ch_assembly)

    ch_reports = METAQUAST_SPADES.out.report_tsv
        .mix(METAQUAST_SPADES.out.report_html)

    ch_logs = SPADES.out.log
        .mix(METAQUAST_SPADES.out.log)

    ch_versions = SPADES.out.versions
        .mix(ch_merge_versions)
        .mix(METAQUAST_SPADES.out.versions)

    emit:
    assembly          = ch_assembly
    individual_assemblies = SPADES.out.contigs
    contig_provenance  = ch_contig_provenance
    assembly_manifest = ch_assembly_manifest
    scaffolds         = SPADES.out.scaffolds
    graph             = SPADES.out.graph
    assembly_log      = SPADES.out.log
    assembly_params   = SPADES.out.params
    metaquast_results = METAQUAST_SPADES.out.results
    metaquast_report  = METAQUAST_SPADES.out.report_tsv
    metaquast_html    = METAQUAST_SPADES.out.report_html
    metaquast_log     = METAQUAST_SPADES.out.log
    reports           = ch_reports
    logs              = ch_logs
    versions          = ch_versions
}
