#!/usr/bin/env nextflow

include { MEGAHIT_ASSEMBLY } from '../megahit_assembly/main'
include { MEGAHIT_BINNING } from '../megahit_binning/main'
include { MEGAHIT_MAG_REFINEMENT } from '../megahit_mag_refinement/main'

workflow MEGAHIT_BRANCH {
    take:
    ch_filtered_reads
    ch_checkm2_db
    ch_gunc_db

    main:
    MEGAHIT_ASSEMBLY(ch_filtered_reads)
    MEGAHIT_BINNING(MEGAHIT_ASSEMBLY.out.assembly, ch_filtered_reads)
    MEGAHIT_MAG_REFINEMENT(
        MEGAHIT_BINNING.out.refined_bins,
        ch_checkm2_db,
        ch_gunc_db
    )

    ch_reports = MEGAHIT_ASSEMBLY.out.reports
        .mix(MEGAHIT_BINNING.out.reports)
        .mix(MEGAHIT_MAG_REFINEMENT.out.reports)

    ch_logs = MEGAHIT_ASSEMBLY.out.logs
        .mix(MEGAHIT_BINNING.out.logs)

    ch_versions = MEGAHIT_ASSEMBLY.out.versions
        .mix(MEGAHIT_BINNING.out.versions)
        .mix(MEGAHIT_MAG_REFINEMENT.out.versions)

    emit:
    assembly            = MEGAHIT_ASSEMBLY.out.assembly
    metaquast_report    = MEGAHIT_ASSEMBLY.out.metaquast_report
    filtered_contigs    = MEGAHIT_BINNING.out.filtered_contigs
    raw_bins            = MEGAHIT_BINNING.out.refined_bins
    dastool_summary     = MEGAHIT_BINNING.out.dastool_summary
    selected_mags       = MEGAHIT_MAG_REFINEMENT.out.selected_mags
    species_reps        = MEGAHIT_MAG_REFINEMENT.out.species_reps
    raw_checkm2         = MEGAHIT_MAG_REFINEMENT.out.raw_checkm2
    raw_gunc            = MEGAHIT_MAG_REFINEMENT.out.raw_gunc
    clean_checkm2       = MEGAHIT_MAG_REFINEMENT.out.clean_checkm2
    clean_gunc          = MEGAHIT_MAG_REFINEMENT.out.clean_gunc
    drep_99_clusters    = MEGAHIT_MAG_REFINEMENT.out.drep_99_clusters
    species_95_clusters = MEGAHIT_MAG_REFINEMENT.out.species_95_clusters
    reports             = ch_reports
    logs                = ch_logs
    versions            = ch_versions
}
