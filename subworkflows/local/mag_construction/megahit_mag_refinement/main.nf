#!/usr/bin/env nextflow

include { MAG_REFINEMENT as MEGAHIT_REFINEMENT_CORE } from '../mag_refinement/main'

workflow MEGAHIT_MAG_REFINEMENT {
    take:
    ch_bins
    ch_checkm2_db
    ch_gunc_db

    main:
    MEGAHIT_REFINEMENT_CORE(ch_bins, ch_checkm2_db, ch_gunc_db)

    emit:
    selected_mags       = MEGAHIT_REFINEMENT_CORE.out.selected_mags
    species_reps        = MEGAHIT_REFINEMENT_CORE.out.species_reps
    raw_checkm2         = MEGAHIT_REFINEMENT_CORE.out.raw_checkm2
    raw_gunc            = MEGAHIT_REFINEMENT_CORE.out.raw_gunc
    clean_checkm2       = MEGAHIT_REFINEMENT_CORE.out.clean_checkm2
    clean_gunc          = MEGAHIT_REFINEMENT_CORE.out.clean_gunc
    hq_table            = MEGAHIT_REFINEMENT_CORE.out.hq_table
    drep_99_clusters    = MEGAHIT_REFINEMENT_CORE.out.drep_99_clusters
    species_95_clusters = MEGAHIT_REFINEMENT_CORE.out.species_95_clusters
    reports             = MEGAHIT_REFINEMENT_CORE.out.reports
    versions            = MEGAHIT_REFINEMENT_CORE.out.versions
}
