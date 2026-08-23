#!/usr/bin/env nextflow

include { CHECKM2 as CHECKM2_RAW } from '../../../../modules/core/checkm2/main'
include { GUNC as GUNC_RAW } from '../../../../modules/core/gunc/main'
include { SELECT_HIGH_QUALITY_MAGS } from '../../../../modules/local/select_high_quality_mags/main'
include { DREP as DREP_ANI_99 } from '../../../../modules/core/drep/main'
include { DREP as DREP_SPECIES_95 } from '../../../../modules/core/drep/main'
include { CHECKM2 as CHECKM2_CLEAN } from '../../../../modules/core/checkm2/main'
include { GUNC as GUNC_CLEAN } from '../../../../modules/core/gunc/main'

workflow MAG_REFINEMENT {
    take:
    ch_bins
    ch_checkm2_db
    ch_gunc_db

    main:
    def pipeline_root = params.pipeline_root ?: projectDir

    CHECKM2_RAW(ch_bins, ch_checkm2_db)
    GUNC_RAW(ch_bins, ch_gunc_db)

    ch_selection_input = ch_bins
        .map { meta, bins -> tuple(meta.id, meta, bins) }
        .join(CHECKM2_RAW.out.quality.map { meta, quality -> tuple(meta.id, quality) })
        .map { _id, meta, bins, quality -> tuple(meta, bins, quality) }

    ch_selection_script = channel.value(
        file("${pipeline_root}/bin/select_high_quality_mags.py", checkIfExists: true)
    )

    SELECT_HIGH_QUALITY_MAGS(
        ch_selection_input,
        ch_selection_script,
        channel.value(params.hq_completeness),
        channel.value(params.hq_contamination)
    )

    DREP_ANI_99(
        SELECT_HIGH_QUALITY_MAGS.out.mags,
        SELECT_HIGH_QUALITY_MAGS.out.genome_info.map { _meta, genome_info -> genome_info },
        channel.value(params.derep_ani),
        channel.value(params.derep_coverage),
        channel.value('ani99')
    )

    DREP_SPECIES_95(
        DREP_ANI_99.out.representatives,
        SELECT_HIGH_QUALITY_MAGS.out.genome_info.map { _meta, genome_info -> genome_info },
        channel.value(params.species_ani),
        channel.value(params.derep_coverage),
        channel.value('species95')
    )

    CHECKM2_CLEAN(DREP_ANI_99.out.representatives, ch_checkm2_db)
    GUNC_CLEAN(DREP_ANI_99.out.representatives, ch_gunc_db)

    ch_reports = CHECKM2_RAW.out.quality
        .mix(GUNC_RAW.out.summary)
        .mix(SELECT_HIGH_QUALITY_MAGS.out.table)
        .mix(DREP_ANI_99.out.clusters)
        .mix(DREP_SPECIES_95.out.clusters)
        .mix(CHECKM2_CLEAN.out.quality)
        .mix(GUNC_CLEAN.out.summary)

    ch_versions = CHECKM2_RAW.out.versions
        .mix(GUNC_RAW.out.versions)
        .mix(SELECT_HIGH_QUALITY_MAGS.out.versions)
        .mix(DREP_ANI_99.out.versions)
        .mix(DREP_SPECIES_95.out.versions)
        .mix(CHECKM2_CLEAN.out.versions)
        .mix(GUNC_CLEAN.out.versions)

    emit:
    selected_mags       = DREP_ANI_99.out.representatives
    species_reps        = DREP_SPECIES_95.out.representatives
    raw_checkm2         = CHECKM2_RAW.out.quality
    raw_gunc            = GUNC_RAW.out.summary
    clean_checkm2       = CHECKM2_CLEAN.out.quality
    clean_gunc          = GUNC_CLEAN.out.summary
    hq_table            = SELECT_HIGH_QUALITY_MAGS.out.table
    drep_99_clusters    = DREP_ANI_99.out.clusters
    species_95_clusters = DREP_SPECIES_95.out.clusters
    reports             = ch_reports
    versions            = ch_versions
}
