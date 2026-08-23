#!/usr/bin/env nextflow

include { MEGAHIT } from '../../../../modules/core/megahit/main'
include { METAQUAST as METAQUAST_MEGAHIT } from '../../../../modules/core/metaquast/main'

workflow MEGAHIT_ASSEMBLY {
    take:
    ch_filtered_reads

    main:
    ch_coassembly_input = ch_filtered_reads
        .collect(flat: false)
        .map { records ->
            if (!records) {
                error "MEGAHIT coassembly requires at least one paired-end sample"
            }

            def ordered_records = records.toList().sort { left, right ->
                left[0]['id'].toString() <=> right[0]['id'].toString()
            }

            ordered_records.each { record ->
                if (record.size() != 2 || !record[0]['id'] || record[1].size() != 2) {
                    error "MEGAHIT coassembly expects tuples of [meta, [read_1, read_2]]"
                }
            }

            def sample_ids = ordered_records.collect { record -> record[0]['id'].toString() }
            def reads      = ordered_records.collectMany { record -> record[1].toList() }
            def meta       = [
                id         : 'megahit_coassembly',
                assembler  : 'megahit',
                branch     : 'megahit',
                sample_ids : sample_ids
            ]

            tuple(meta, reads)
        }

    MEGAHIT(ch_coassembly_input)
    METAQUAST_MEGAHIT(MEGAHIT.out.contigs)

    ch_reports = METAQUAST_MEGAHIT.out.report_tsv
        .mix(METAQUAST_MEGAHIT.out.report_html)

    ch_logs = MEGAHIT.out.log
        .mix(METAQUAST_MEGAHIT.out.log)

    ch_versions = MEGAHIT.out.versions
        .mix(METAQUAST_MEGAHIT.out.versions)

    emit:
    assembly          = MEGAHIT.out.contigs
    assembly_log      = MEGAHIT.out.log
    metaquast_results = METAQUAST_MEGAHIT.out.results
    metaquast_report  = METAQUAST_MEGAHIT.out.report_tsv
    metaquast_html    = METAQUAST_MEGAHIT.out.report_html
    metaquast_log     = METAQUAST_MEGAHIT.out.log
    reports           = ch_reports
    logs              = ch_logs
    versions          = ch_versions
}
