#!/usr/bin/env nextflow

include { MEGAHIT_BINNING } from '../../subworkflows/local/mag_construction/megahit_binning/main'
include { SPADES_BINNING } from '../../subworkflows/local/mag_construction/spades_binning/main'

workflow {
    def pipeline_root = "${projectDir}/../.."
    def generated = "${pipeline_root}/tests/generated_data"

    ch_reads = channel.of(
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

    ch_megahit_assembly = channel.value(
        tuple(
            [id: 'megahit_coassembly', assembler: 'megahit', branch: 'megahit', sample_ids: ['sample_A', 'sample_B']],
            file("${generated}/microbial_fixture.fa", checkIfExists: true)
        )
    )
    ch_spades_assembly = channel.value(
        tuple(
            [id: 'spades_coassembly', assembler: 'spades', branch: 'spades', sample_ids: ['sample_A', 'sample_B']],
            file("${generated}/microbial_fixture.fa", checkIfExists: true)
        )
    )

    MEGAHIT_BINNING(ch_megahit_assembly, ch_reads)
    SPADES_BINNING(ch_spades_assembly, ch_reads)
}
