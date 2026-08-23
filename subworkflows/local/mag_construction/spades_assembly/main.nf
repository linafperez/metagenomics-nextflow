#!/usr/bin/env nextflow

include { SPADES } from '../../../../modules/core/spades/main'
include { METAQUAST as METAQUAST_SPADES } from '../../../../modules/core/metaquast/main'

workflow SPADES_ASSEMBLY {
    take:
    ch_filtered_reads

    main:
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

    SPADES(ch_coassembly_input)
    METAQUAST_SPADES(SPADES.out.contigs)

    ch_reports = METAQUAST_SPADES.out.report_tsv
        .mix(METAQUAST_SPADES.out.report_html)

    ch_logs = SPADES.out.log
        .mix(METAQUAST_SPADES.out.log)

    ch_versions = SPADES.out.versions
        .mix(METAQUAST_SPADES.out.versions)

    emit:
    assembly          = SPADES.out.contigs
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
