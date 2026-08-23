#!/usr/bin/env nextflow

include { MAG_CONSTRUCTION } from '../../subworkflows/local/mag_construction/main'

workflow {
    def pipeline_root = "${projectDir}/../.."
    def generated = "${pipeline_root}/tests/generated_data"

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

    ch_checkm2_db = channel.value(
        file("${generated}/databases/checkm2/uniref100.KO.1.dmnd", checkIfExists: true)
    )
    ch_gunc_db = channel.value(
        file("${generated}/databases/gunc/gunc_db.dmnd", checkIfExists: true)
    )

    MAG_CONSTRUCTION(ch_filtered_reads, ch_checkm2_db, ch_gunc_db)
}
