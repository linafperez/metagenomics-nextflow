#!/usr/bin/env nextflow

include { CHECK_SAMPLESHEET } from '../../modules/local/check_samplesheet/main'
include { MEGAHIT_ASSEMBLY } from '../../subworkflows/local/mag_construction/megahit_assembly/main'
include { SPADES_ASSEMBLY } from '../../subworkflows/local/mag_construction/spades_assembly/main'
include { BENCHMARK_BINNING as MEGAHIT_BENCHMARK_BINNING } from '../subworkflows/benchmark_binning/main'
include { BENCHMARK_BINNING as SPADES_BENCHMARK_BINNING } from '../subworkflows/benchmark_binning/main'
include { MEGAHIT_MAG_REFINEMENT } from '../../subworkflows/local/mag_construction/megahit_mag_refinement/main'
include { SPADES_MAG_REFINEMENT } from '../../subworkflows/local/mag_construction/spades_mag_refinement/main'
include { FINAL_MAG_CATALOG } from '../../subworkflows/local/mag_construction/final_catalog/main'

workflow BENCHMARK {
    main:
    def assembler = (params.benchmark_assembler ?: '').toString().toLowerCase()
    def binner = (params.benchmark_binner ?: '').toString().toLowerCase()
    def valid_assemblers = ['megahit', 'spades', 'both']
    def valid_binners = ['comebin', 'metabat2', 'semibin2', 'vamb', 'all']

    if (!valid_assemblers.contains(assembler)) {
        error "Unsupported --benchmark_assembler '${params.benchmark_assembler}'. Expected one of: ${valid_assemblers.join(', ')}"
    }
    if (!valid_binners.contains(binner)) {
        error "Unsupported --benchmark_binner '${params.benchmark_binner}'. Expected one of: ${valid_binners.join(', ')}"
    }
    if (!params.input) {
        error 'Missing required parameter: --input'
    }
    if (!params.checkm2_db) {
        error 'Missing required parameter: --checkm2_db'
    }
    if (!params.gunc_db) {
        error 'Missing required parameter: --gunc_db'
    }

    def pipeline_root = params.pipeline_root ?: "${projectDir}/../.."
    ch_samplesheet = channel.fromPath(params.input, checkIfExists: true)
    ch_validator = channel.value(
        file("${pipeline_root}/bin/check_samplesheet.py", checkIfExists: true)
    )
    CHECK_SAMPLESHEET(ch_samplesheet, ch_validator)

    ch_filtered_reads = CHECK_SAMPLESHEET.out.csv
        .splitCsv(header: true)
        .map { row ->
            tuple(
                [id: row.sample, single_end: false],
                [
                    file(row.fastq_1, checkIfExists: true),
                    file(row.fastq_2, checkIfExists: true)
                ]
            )
        }

    ch_checkm2_db = channel.value(file(params.checkm2_db, checkIfExists: true))
    ch_gunc_db = channel.value(file(params.gunc_db, checkIfExists: true))

    if (assembler == 'megahit' || assembler == 'both') {
        MEGAHIT_ASSEMBLY(ch_filtered_reads)
        MEGAHIT_BENCHMARK_BINNING(
            MEGAHIT_ASSEMBLY.out.assembly,
            ch_filtered_reads,
            binner
        )
        MEGAHIT_MAG_REFINEMENT(
            MEGAHIT_BENCHMARK_BINNING.out.bins,
            ch_checkm2_db,
            ch_gunc_db
        )
    }

    if (assembler == 'spades' || assembler == 'both') {
        SPADES_ASSEMBLY(ch_filtered_reads)
        SPADES_BENCHMARK_BINNING(
            SPADES_ASSEMBLY.out.assembly,
            ch_filtered_reads,
            binner
        )
        SPADES_MAG_REFINEMENT(
            SPADES_BENCHMARK_BINNING.out.bins,
            ch_checkm2_db,
            ch_gunc_db
        )
    }

    if (assembler == 'both') {
        FINAL_MAG_CATALOG(
            MEGAHIT_MAG_REFINEMENT.out.selected_mags,
            MEGAHIT_MAG_REFINEMENT.out.clean_checkm2,
            SPADES_MAG_REFINEMENT.out.selected_mags,
            SPADES_MAG_REFINEMENT.out.clean_checkm2,
            ch_checkm2_db,
            ch_gunc_db
        )
    }
}

workflow {
    BENCHMARK()
}
