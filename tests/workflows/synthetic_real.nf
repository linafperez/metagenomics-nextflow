#!/usr/bin/env nextflow

include { CHECK_SAMPLESHEET } from '../../modules/local/check_samplesheet/main'
include { BOWTIE2_BUILD_SYNTHETIC } from '../modules/bowtie2_build/main'
include { QUALITY_CONTROL_AND_FILTERING } from '../../subworkflows/local/quality_control_and_filtering/main'
include { MEGAHIT_ASSEMBLY } from '../../subworkflows/local/mag_construction/megahit_assembly/main'
include { SPADES_ASSEMBLY } from '../../subworkflows/local/mag_construction/spades_assembly/main'

workflow SYNTHETIC_REAL_TOOLS {
    main:
    def pipeline_root = params.pipeline_root ?: projectDir

    if (params.input == null || params.input.toString().trim().isEmpty()) {
        error 'Synthetic real-tool validation requires --input'
    }

    ch_samplesheet = channel.fromPath(params.input, checkIfExists: true)
    ch_samplesheet_validator = channel.value(
        file("${pipeline_root}/bin/check_samplesheet.py", checkIfExists: true)
    )
    ch_host_reference = channel.value(
        file("${pipeline_root}/tests/generated_reference/GRCh38.p14.fa", checkIfExists: true)
    )

    CHECK_SAMPLESHEET(ch_samplesheet, ch_samplesheet_validator)
    BOWTIE2_BUILD_SYNTHETIC(ch_host_reference)

    CHECK_SAMPLESHEET.out.csv
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample, single_end: false]
            def reads = [
                file(row.fastq_1, checkIfExists: true),
                file(row.fastq_2, checkIfExists: true)
            ]
            tuple(meta, reads)
        }
        .set { ch_raw_reads }

    QUALITY_CONTROL_AND_FILTERING(
        ch_raw_reads,
        BOWTIE2_BUILD_SYNTHETIC.out.index,
        BOWTIE2_BUILD_SYNTHETIC.out.prefix
    )

    MEGAHIT_ASSEMBLY(QUALITY_CONTROL_AND_FILTERING.out.filtered_reads)
    SPADES_ASSEMBLY(QUALITY_CONTROL_AND_FILTERING.out.filtered_reads)

    ch_versions = CHECK_SAMPLESHEET.out.versions
        .mix(BOWTIE2_BUILD_SYNTHETIC.out.versions)
        .mix(QUALITY_CONTROL_AND_FILTERING.out.versions)
        .mix(MEGAHIT_ASSEMBLY.out.versions)
        .mix(SPADES_ASSEMBLY.out.versions)

    emit:
    filtered_reads     = QUALITY_CONTROL_AND_FILTERING.out.filtered_reads
    megahit_assembly   = MEGAHIT_ASSEMBLY.out.assembly
    spades_assembly    = SPADES_ASSEMBLY.out.assembly
    megahit_metaquast  = MEGAHIT_ASSEMBLY.out.metaquast_report
    spades_metaquast   = SPADES_ASSEMBLY.out.metaquast_report
    versions           = ch_versions
}

workflow {
    SYNTHETIC_REAL_TOOLS()
}
