#!/usr/bin/env nextflow

include { COMBINE_MAG_CATALOGS } from '../../../../modules/local/combine_mag_catalogs/main'
include { DREP as DREP_FINAL_99 } from '../../../../modules/core/drep/main'
include { DREP as DREP_FINAL_SPECIES_95 } from '../../../../modules/core/drep/main'
include { FINALIZE_MAG_CATALOG } from '../../../../modules/local/finalize_mag_catalog/main'
include { CHECKM2 as CHECKM2_FINAL } from '../../../../modules/core/checkm2/main'
include { GUNC as GUNC_FINAL } from '../../../../modules/core/gunc/main'

workflow FINAL_MAG_CATALOG {
    take:
    ch_megahit_mags
    ch_megahit_quality
    ch_spades_mags
    ch_spades_quality
    ch_checkm2_db
    ch_gunc_db

    main:
    ch_megahit_keyed = ch_megahit_mags
        .map { meta, mags -> tuple(meta.id, meta, mags) }
        .join(ch_megahit_quality.map { meta, quality -> tuple(meta.id, quality) })

    ch_spades_keyed = ch_spades_mags
        .map { meta, mags -> tuple(meta.id, meta, mags) }
        .join(ch_spades_quality.map { meta, quality -> tuple(meta.id, quality) })

    ch_combined_input = ch_megahit_keyed
        .combine(ch_spades_keyed)
        .map { _megahit_id, _megahit_meta, megahit_mags, megahit_quality, _spades_id, _spades_meta, spades_mags, spades_quality ->
            def meta = [id: 'final_catalog', branch: 'final', assembler: 'combined']
            tuple(meta, megahit_mags, megahit_quality, spades_mags, spades_quality)
        }

    ch_combine_script = channel.value(
        file("${projectDir}/bin/combine_mag_catalogs.py", checkIfExists: true)
    )
    ch_finalize_script = channel.value(
        file("${projectDir}/bin/finalize_mag_catalog.py", checkIfExists: true)
    )

    COMBINE_MAG_CATALOGS(ch_combined_input, ch_combine_script)

    DREP_FINAL_99(
        COMBINE_MAG_CATALOGS.out.mags,
        COMBINE_MAG_CATALOGS.out.genome_info.map { _meta, genome_info -> genome_info },
        channel.value(params.derep_ani),
        channel.value(params.derep_coverage),
        channel.value('ani99')
    )

    DREP_FINAL_SPECIES_95(
        DREP_FINAL_99.out.representatives,
        COMBINE_MAG_CATALOGS.out.genome_info.map { _meta, genome_info -> genome_info },
        channel.value(params.species_ani),
        channel.value(params.derep_coverage),
        channel.value('species95')
    )

    ch_finalize_input = DREP_FINAL_99.out.representatives
        .map { meta, mags -> tuple(meta.id, meta, mags) }
        .join(COMBINE_MAG_CATALOGS.out.provenance.map { meta, provenance -> tuple(meta.id, provenance) })
        .join(COMBINE_MAG_CATALOGS.out.quality.map { meta, quality -> tuple(meta.id, quality) })
        .map { _id, meta, mags, provenance, quality -> tuple(meta, mags, provenance, quality) }

    FINALIZE_MAG_CATALOG(ch_finalize_input, ch_finalize_script)

    CHECKM2_FINAL(FINALIZE_MAG_CATALOG.out.mags, ch_checkm2_db)
    GUNC_FINAL(FINALIZE_MAG_CATALOG.out.mags, ch_gunc_db)

    ch_reports = CHECKM2_FINAL.out.quality
        .mix(GUNC_FINAL.out.summary)
        .mix(COMBINE_MAG_CATALOGS.out.quality)
        .mix(COMBINE_MAG_CATALOGS.out.provenance)
        .mix(FINALIZE_MAG_CATALOG.out.quality)
        .mix(FINALIZE_MAG_CATALOG.out.provenance)
        .mix(DREP_FINAL_99.out.clusters)
        .mix(DREP_FINAL_SPECIES_95.out.clusters)

    ch_versions = COMBINE_MAG_CATALOGS.out.versions
        .mix(DREP_FINAL_99.out.versions)
        .mix(DREP_FINAL_SPECIES_95.out.versions)
        .mix(FINALIZE_MAG_CATALOG.out.versions)
        .mix(CHECKM2_FINAL.out.versions)
        .mix(GUNC_FINAL.out.versions)

    emit:
    mags                = FINALIZE_MAG_CATALOG.out.mags
    provenance          = FINALIZE_MAG_CATALOG.out.provenance
    quality             = FINALIZE_MAG_CATALOG.out.quality
    species_reps        = DREP_FINAL_SPECIES_95.out.representatives
    drep_99_clusters    = DREP_FINAL_99.out.clusters
    species_95_clusters = DREP_FINAL_SPECIES_95.out.clusters
    checkm2             = CHECKM2_FINAL.out.quality
    gunc                = GUNC_FINAL.out.summary
    reports             = ch_reports
    versions            = ch_versions
}
