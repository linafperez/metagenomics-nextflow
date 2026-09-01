#!/usr/bin/env nextflow

include { CHECK_SAMPLESHEET } from '../modules/local/check_samplesheet/main'

include { QUALITY_CONTROL_AND_FILTERING } from '../subworkflows/local/quality_control_and_filtering/main'
include { METAGENOMICS_GLOBAL } from '../subworkflows/local/metagenomics_global/main'

workflow METAGENOMICS {
    main:
    if (params.sraProject) {
        error '--input and --sra-project are mutually exclusive; use the launcher for SRA project mode'
    }

    def required_parameters = [
        input              : params.input,
        host_bowtie2_index : params.host_bowtie2_index
    ]

    def missing_parameters = required_parameters
        .findAll { _name, value -> value == null || value.toString().trim().isEmpty() }
        .keySet()
        .sort()

    if (missing_parameters) {
        error "Missing required parameter(s): ${missing_parameters.collect { name -> "--${name}" }.join(', ')}"
    }

    ch_samplesheet = channel.fromPath(params.input, checkIfExists: true)
    ch_samplesheet_validator = channel.value(
        file("${projectDir}/bin/check_samplesheet.py", checkIfExists: true)
    )
    CHECK_SAMPLESHEET(ch_samplesheet, ch_samplesheet_validator)

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

    def host_index_prefix = file(params.host_bowtie2_index).name
    ch_host_index = channel
        .fromPath("${params.host_bowtie2_index}*.bt2*", checkIfExists: true)
        .collect()
        .map { index_files ->
            if (index_files.size() != 6) {
                error "Bowtie2 index prefix '${params.host_bowtie2_index}' resolved to ${index_files.size()} files; exactly six are required"
            }
            index_files.toList().sort { left, right -> left.name <=> right.name }
        }

    ch_host_index_prefix = channel.value(host_index_prefix)
    QUALITY_CONTROL_AND_FILTERING(
        ch_raw_reads,
        ch_host_index,
        ch_host_index_prefix
    )

    ch_preprocessing_reports = QUALITY_CONTROL_AND_FILTERING.out.raw_fastqc
        .mix(QUALITY_CONTROL_AND_FILTERING.out.fastp_json)
        .mix(QUALITY_CONTROL_AND_FILTERING.out.clean_fastqc)
        .mix(QUALITY_CONTROL_AND_FILTERING.out.bowtie2_logs)

    ch_preprocessing_versions = CHECK_SAMPLESHEET.out.versions
        .mix(QUALITY_CONTROL_AND_FILTERING.out.versions)

    METAGENOMICS_GLOBAL(
        QUALITY_CONTROL_AND_FILTERING.out.filtered_reads,
        ch_preprocessing_reports,
        ch_preprocessing_versions
    )

    emit:
    filtered_reads                  = QUALITY_CONTROL_AND_FILTERING.out.filtered_reads
    final_mags                      = METAGENOMICS_GLOBAL.out.final_mags
    final_catalog_provenance        = METAGENOMICS_GLOBAL.out.final_catalog_provenance
    final_catalog_quality           = METAGENOMICS_GLOBAL.out.final_catalog_quality
    final_checkm2                   = METAGENOMICS_GLOBAL.out.final_checkm2
    final_gunc                      = METAGENOMICS_GLOBAL.out.final_gunc
    final_species_representatives   = METAGENOMICS_GLOBAL.out.final_species_representatives
    gtdbtk_summaries                = METAGENOMICS_GLOBAL.out.gtdbtk_summaries
    phylogenomic_tree               = METAGENOMICS_GLOBAL.out.phylogenomic_tree
    functional_annotations         = METAGENOMICS_GLOBAL.out.functional_annotations
    mag_abundance                   = METAGENOMICS_GLOBAL.out.mag_abundance
    multiqc_report                  = METAGENOMICS_GLOBAL.out.multiqc_report
    software_versions               = METAGENOMICS_GLOBAL.out.software_versions
}
