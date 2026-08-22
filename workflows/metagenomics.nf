#!/usr/bin/env nextflow

include { CHECK_SAMPLESHEET } from '../modules/local/check_samplesheet/main'

include { QUALITY_CONTROL_AND_FILTERING } from '../subworkflows/local/quality_control_and_filtering/main'
include { MAG_CONSTRUCTION } from '../subworkflows/local/mag_construction/main'
include { TAXONOMIC_CLASSIFICATION_AND_PHYLOGENOMICS } from '../subworkflows/local/taxonomic_classification_and_phylogenomics/main'
include { GENE_PREDICTION_AND_FUNCTIONAL_ANNOTATION } from '../subworkflows/local/gene_prediction_and_functional_annotation/main'
include { MAG_ABUNDANCE_ESTIMATION } from '../subworkflows/local/mag_abundance_estimation/main'
include { GLOBAL_PROCESSING_EVALUATION } from '../subworkflows/local/global_processing_evaluation/main'

workflow METAGENOMICS {
    main:
    if (!params.input) {
        error "Missing required parameter: --input"
    }

    if (!params.host_bowtie2_index) {
        error "Missing required parameter: --host_bowtie2_index"
    }

    ch_samplesheet = channel.fromPath(params.input, checkIfExists: true)
    ch_samplesheet_validator = channel.value(
        file("${projectDir}/bin/check_samplesheet.py", checkIfExists: true)
    )
    CHECK_SAMPLESHEET(ch_samplesheet, ch_samplesheet_validator)

    CHECK_SAMPLESHEET.out.csv
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                id         : row.sample,
                single_end : false
            ]

            def reads = [
                file(row.fastq_1, checkIfExists: true),
                file(row.fastq_2, checkIfExists: true)
            ]

            tuple(meta, reads)
        }
        .set { ch_raw_reads }

    def host_index_prefix = file(params.host_bowtie2_index).name

    ch_host_index = channel
        .fromPath("${params.host_bowtie2_index}*.bt2*", checkIfExists: true)
        .collect()
        .map { index_files ->
            if (index_files.size() != 6) {
                error "Bowtie2 index prefix '${params.host_bowtie2_index}' resolved to ${index_files.size()} files; exactly 6 are required"
            }
            index_files
        }

    ch_host_index_prefix = channel.value(host_index_prefix)

    QUALITY_CONTROL_AND_FILTERING(
        ch_raw_reads,
        ch_host_index,
        ch_host_index_prefix
    )

    // Phase 1 intentionally stops here. The remaining scientific subworkflows
    // are represented as skeletons but are not invoked until later phases.

    emit:
    filtered_reads = QUALITY_CONTROL_AND_FILTERING.out.filtered_reads
    raw_fastqc     = QUALITY_CONTROL_AND_FILTERING.out.raw_fastqc
    fastp_json     = QUALITY_CONTROL_AND_FILTERING.out.fastp_json
    clean_fastqc   = QUALITY_CONTROL_AND_FILTERING.out.clean_fastqc
    bowtie2_logs   = QUALITY_CONTROL_AND_FILTERING.out.bowtie2_logs
    versions       = QUALITY_CONTROL_AND_FILTERING.out.versions.mix(CHECK_SAMPLESHEET.out.versions)
}
