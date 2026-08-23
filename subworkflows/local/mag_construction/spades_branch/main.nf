#!/usr/bin/env nextflow

include { SPADES_ASSEMBLY } from '../spades_assembly/main'
include { SPADES_BINNING } from '../spades_binning/main'
include { SPADES_MAG_REFINEMENT } from '../spades_mag_refinement/main'

workflow SPADES_BRANCH {
    take:
    ch_filtered_reads
    ch_checkm2_db
    ch_gunc_db

    main:
    SPADES_ASSEMBLY(ch_filtered_reads)
    SPADES_BINNING(SPADES_ASSEMBLY.out.assembly, ch_filtered_reads)
    SPADES_MAG_REFINEMENT(
        SPADES_BINNING.out.refined_bins,
        ch_checkm2_db,
        ch_gunc_db
    )

    ch_reports = SPADES_ASSEMBLY.out.reports
        .mix(SPADES_BINNING.out.reports)
        .mix(SPADES_MAG_REFINEMENT.out.reports)

    ch_logs = SPADES_ASSEMBLY.out.logs
        .mix(SPADES_BINNING.out.logs)

    ch_versions = SPADES_ASSEMBLY.out.versions
        .mix(SPADES_BINNING.out.versions)
        .mix(SPADES_MAG_REFINEMENT.out.versions)

    emit:
    assembly            = SPADES_ASSEMBLY.out.assembly
    metaquast_report    = SPADES_ASSEMBLY.out.metaquast_report
    filtered_contigs    = SPADES_BINNING.out.filtered_contigs
    raw_bins            = SPADES_BINNING.out.refined_bins
    dastool_summary     = SPADES_BINNING.out.dastool_summary
    selected_mags       = SPADES_MAG_REFINEMENT.out.selected_mags
    species_reps        = SPADES_MAG_REFINEMENT.out.species_reps
    raw_checkm2         = SPADES_MAG_REFINEMENT.out.raw_checkm2
    raw_gunc            = SPADES_MAG_REFINEMENT.out.raw_gunc
    clean_checkm2       = SPADES_MAG_REFINEMENT.out.clean_checkm2
    clean_gunc          = SPADES_MAG_REFINEMENT.out.clean_gunc
    drep_99_clusters    = SPADES_MAG_REFINEMENT.out.drep_99_clusters
    species_95_clusters = SPADES_MAG_REFINEMENT.out.species_95_clusters
    reports             = ch_reports
    logs                = ch_logs
    versions            = ch_versions
}
