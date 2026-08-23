#!/usr/bin/env nextflow

include { MEGAHIT_BRANCH } from './megahit_branch/main'
include { SPADES_BRANCH } from './spades_branch/main'
include { FINAL_MAG_CATALOG } from './final_catalog/main'

workflow MAG_CONSTRUCTION {
    take:
    ch_filtered_reads
    ch_checkm2_db
    ch_gunc_db

    main:
    MEGAHIT_BRANCH(ch_filtered_reads, ch_checkm2_db, ch_gunc_db)
    SPADES_BRANCH(ch_filtered_reads, ch_checkm2_db, ch_gunc_db)

    FINAL_MAG_CATALOG(
        MEGAHIT_BRANCH.out.selected_mags,
        MEGAHIT_BRANCH.out.clean_checkm2,
        SPADES_BRANCH.out.selected_mags,
        SPADES_BRANCH.out.clean_checkm2,
        ch_checkm2_db,
        ch_gunc_db
    )

    ch_reports = MEGAHIT_BRANCH.out.reports
        .mix(SPADES_BRANCH.out.reports)
        .mix(FINAL_MAG_CATALOG.out.reports)

    ch_logs = MEGAHIT_BRANCH.out.logs
        .mix(SPADES_BRANCH.out.logs)

    ch_versions = MEGAHIT_BRANCH.out.versions
        .mix(SPADES_BRANCH.out.versions)
        .mix(FINAL_MAG_CATALOG.out.versions)

    emit:
    final_mags          = FINAL_MAG_CATALOG.out.mags
    final_provenance    = FINAL_MAG_CATALOG.out.provenance
    final_quality       = FINAL_MAG_CATALOG.out.quality
    final_checkm2       = FINAL_MAG_CATALOG.out.checkm2
    final_gunc          = FINAL_MAG_CATALOG.out.gunc
    final_species_reps  = FINAL_MAG_CATALOG.out.species_reps
    final_99_clusters   = FINAL_MAG_CATALOG.out.drep_99_clusters
    final_95_clusters   = FINAL_MAG_CATALOG.out.species_95_clusters
    megahit_assembly    = MEGAHIT_BRANCH.out.assembly
    megahit_metaquast   = MEGAHIT_BRANCH.out.metaquast_report
    megahit_raw_bins    = MEGAHIT_BRANCH.out.raw_bins
    megahit_clean_mags  = MEGAHIT_BRANCH.out.selected_mags
    spades_assembly     = SPADES_BRANCH.out.assembly
    spades_metaquast    = SPADES_BRANCH.out.metaquast_report
    spades_raw_bins     = SPADES_BRANCH.out.raw_bins
    spades_clean_mags   = SPADES_BRANCH.out.selected_mags
    reports             = ch_reports
    logs                = ch_logs
    versions            = ch_versions
}
