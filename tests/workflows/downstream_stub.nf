#!/usr/bin/env nextflow

include { GENE_PREDICTION_AND_FUNCTIONAL_ANNOTATION } from '../../subworkflows/local/gene_prediction_and_functional_annotation/main'
include { MAG_ABUNDANCE_ESTIMATION } from '../../subworkflows/local/mag_abundance_estimation/main'

workflow {
    def pipeline_root = "${projectDir}/../.."
    def generated = "${pipeline_root}/tests/generated_data"

    ch_final_mags = channel.value(
        tuple(
            [id: 'final_catalog', catalog_id: 'final_catalog', branch: 'combined'],
            [
                file("${generated}/bins/megahit/bin_001.fa", checkIfExists: true),
                file("${generated}/bins/spades/bin_002.fa", checkIfExists: true)
            ]
        )
    )

    ch_filtered_reads = channel.of(
        tuple(
            [id: 'sample_A', single_end: false],
            [
                file("${generated}/sample_A_R1.fastq.gz", checkIfExists: true),
                file("${generated}/sample_A_R2.fastq.gz", checkIfExists: true)
            ]
        ),
        tuple(
            [id: 'sample_B', single_end: false],
            [
                file("${generated}/sample_B_R1.fastq.gz", checkIfExists: true),
                file("${generated}/sample_B_R2.fastq.gz", checkIfExists: true)
            ]
        )
    )

    ch_genemark_home = channel.value(file("${generated}/software/genemark", checkIfExists: true))
    ch_genemark_key = channel.value(file("${generated}/licenses/gm_key", checkIfExists: true))
    ch_eggnog_db = channel.value(file("${generated}/databases/eggnog", checkIfExists: true))
    ch_interproscan_data = channel.value(file("${generated}/databases/interproscan", checkIfExists: true))

    GENE_PREDICTION_AND_FUNCTIONAL_ANNOTATION(
        ch_final_mags,
        ch_genemark_home,
        ch_genemark_key,
        ch_eggnog_db,
        ch_interproscan_data
    )

    MAG_ABUNDANCE_ESTIMATION(ch_final_mags, ch_filtered_reads)
}
