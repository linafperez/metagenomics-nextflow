#!/usr/bin/env nextflow

include { MEGAHIT_MAG_REFINEMENT } from '../../subworkflows/local/mag_construction/megahit_mag_refinement/main'
include { SPADES_MAG_REFINEMENT } from '../../subworkflows/local/mag_construction/spades_mag_refinement/main'
include { FINAL_MAG_CATALOG } from '../../subworkflows/local/mag_construction/final_catalog/main'

workflow {
    def pipeline_root = "${projectDir}/../.."
    def generated = "${pipeline_root}/tests/generated_data"

    ch_megahit_bins = channel.value(
        tuple(
            [id: 'megahit_bins', assembler: 'megahit', branch: 'megahit'],
            [
                file("${generated}/bins/megahit/bin_001.fa", checkIfExists: true),
                file("${generated}/bins/megahit/bin_002.fa", checkIfExists: true)
            ]
        )
    )

    ch_spades_bins = channel.value(
        tuple(
            [id: 'spades_bins', assembler: 'spades', branch: 'spades'],
            [
                file("${generated}/bins/spades/bin_001.fa", checkIfExists: true),
                file("${generated}/bins/spades/bin_002.fa", checkIfExists: true)
            ]
        )
    )

    ch_checkm2_db = channel.value(
        file("${generated}/databases/checkm2/uniref100.KO.1.dmnd", checkIfExists: true)
    )
    ch_gunc_db = channel.value(
        file("${generated}/databases/gunc/gunc_db.dmnd", checkIfExists: true)
    )

    MEGAHIT_MAG_REFINEMENT(ch_megahit_bins, ch_checkm2_db, ch_gunc_db)
    SPADES_MAG_REFINEMENT(ch_spades_bins, ch_checkm2_db, ch_gunc_db)

    FINAL_MAG_CATALOG(
        MEGAHIT_MAG_REFINEMENT.out.selected_mags,
        MEGAHIT_MAG_REFINEMENT.out.clean_checkm2,
        SPADES_MAG_REFINEMENT.out.selected_mags,
        SPADES_MAG_REFINEMENT.out.clean_checkm2,
        ch_checkm2_db,
        ch_gunc_db
    )
}
