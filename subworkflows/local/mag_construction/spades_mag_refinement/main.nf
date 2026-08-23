#!/usr/bin/env nextflow

include { MAG_REFINEMENT as SPADES_REFINEMENT_CORE } from '../mag_refinement/main'

workflow SPADES_MAG_REFINEMENT {
    take:
    ch_bins
    ch_checkm2_db
    ch_gunc_db

    main:
    SPADES_REFINEMENT_CORE(ch_bins, ch_checkm2_db, ch_gunc_db)

    emit:
    selected_mags       = SPADES_REFINEMENT_CORE.out.selected_mags
    species_reps        = SPADES_REFINEMENT_CORE.out.species_reps
    raw_checkm2         = SPADES_REFINEMENT_CORE.out.raw_checkm2
    raw_gunc            = SPADES_REFINEMENT_CORE.out.raw_gunc
    clean_checkm2       = SPADES_REFINEMENT_CORE.out.clean_checkm2
    clean_gunc          = SPADES_REFINEMENT_CORE.out.clean_gunc
    hq_table            = SPADES_REFINEMENT_CORE.out.hq_table
    drep_99_clusters    = SPADES_REFINEMENT_CORE.out.drep_99_clusters
    species_95_clusters = SPADES_REFINEMENT_CORE.out.species_95_clusters
    reports             = SPADES_REFINEMENT_CORE.out.reports
    versions            = SPADES_REFINEMENT_CORE.out.versions
}
