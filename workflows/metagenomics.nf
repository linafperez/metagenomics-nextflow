#!/usr/bin/env nextflow

include { CHECK_SAMPLESHEET } from '../modules/local/check_samplesheet/main'
include { COLLECT_VERSIONS } from '../modules/local/collect_versions/main'

include { QUALITY_CONTROL_AND_FILTERING } from '../subworkflows/local/quality_control_and_filtering/main'
include { MAG_CONSTRUCTION } from '../subworkflows/local/mag_construction/main'
include { TAXONOMIC_CLASSIFICATION_AND_PHYLOGENOMICS } from '../subworkflows/local/taxonomic_classification_and_phylogenomics/main'
include { GENE_PREDICTION_AND_FUNCTIONAL_ANNOTATION } from '../subworkflows/local/gene_prediction_and_functional_annotation/main'
include { MAG_ABUNDANCE_ESTIMATION } from '../subworkflows/local/mag_abundance_estimation/main'
include { GLOBAL_PROCESSING_EVALUATION } from '../subworkflows/local/global_processing_evaluation/main'

workflow METAGENOMICS {
    main:
    def required_parameters = [
        input              : params.input,
        host_bowtie2_index : params.host_bowtie2_index,
        checkm2_db         : params.checkm2_db,
        gunc_db            : params.gunc_db,
        gtdbtk_db          : params.gtdbtk_db,
        phylophlan_db      : params.phylophlan_db,
        phylophlan_config  : params.phylophlan_config,
        genemark_home      : params.genemark_home,
        genemark_key       : params.genemark_key,
        eggnog_db          : params.eggnog_db,
        interproscan_data  : params.interproscan_data
    ]

    def missing_parameters = required_parameters
        .findAll { _name, value -> value == null || value.toString().trim().isEmpty() }
        .keySet()
        .sort()

    if (missing_parameters) {
        error "Missing required parameter(s): ${missing_parameters.collect { name -> "--${name}" }.join(', ')}"
    }

    def integer_pattern = ~/^[0-9]+$/
    if (!(params.min_contig_length.toString() ==~ integer_pattern) || params.min_contig_length.toInteger() < 1) {
        error '--min_contig_length must be a positive integer'
    }

    def number_pattern = ~/^(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)$/
    def numeric_parameters = [
        hq_completeness  : params.hq_completeness,
        hq_contamination : params.hq_contamination,
        derep_ani        : params.derep_ani,
        species_ani      : params.species_ani,
        derep_coverage   : params.derep_coverage
    ]
    def invalid_numeric_parameters = numeric_parameters.findAll { _name, value ->
        value == null || !(value.toString() ==~ number_pattern)
    }.keySet().sort()
    if (invalid_numeric_parameters) {
        error "Numeric parameter(s) have invalid values: ${invalid_numeric_parameters.collect { name -> "--${name}" }.join(', ')}"
    }

    def hq_completeness = params.hq_completeness.toDouble()
    def hq_contamination = params.hq_contamination.toDouble()
    def derep_ani = params.derep_ani.toDouble()
    def species_ani = params.species_ani.toDouble()
    def derep_coverage = params.derep_coverage.toDouble()

    if (hq_completeness < 0 || hq_completeness > 100) {
        error '--hq_completeness must be between 0 and 100'
    }
    if (hq_contamination < 0 || hq_contamination > 100) {
        error '--hq_contamination must be between 0 and 100'
    }
    if (derep_ani <= 0 || derep_ani > 1 || species_ani <= 0 || species_ani > 1) {
        error '--derep_ani and --species_ani must be greater than 0 and at most 1'
    }
    if (derep_ani <= species_ani) {
        error '--derep_ani must be greater than --species_ani'
    }
    if (derep_coverage <= 0 || derep_coverage > 1) {
        error '--derep_coverage must be greater than 0 and at most 1'
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
    ch_checkm2_db = channel.value(file(params.checkm2_db, checkIfExists: true))
    ch_gunc_db = channel.value(file(params.gunc_db, checkIfExists: true))
    ch_gtdbtk_db = channel.value(file(params.gtdbtk_db, checkIfExists: true))
    ch_phylophlan_db = channel.value(file(params.phylophlan_db, checkIfExists: true))
    ch_phylophlan_config = channel.value(file(params.phylophlan_config, checkIfExists: true))
    ch_genemark_home = channel.value(file(params.genemark_home, checkIfExists: true))
    ch_genemark_key = channel.value(file(params.genemark_key, checkIfExists: true))
    ch_eggnog_db = channel.value(file(params.eggnog_db, checkIfExists: true))
    ch_interproscan_data = channel.value(file(params.interproscan_data, checkIfExists: true))
    ch_multiqc_config = channel.value(
        file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)
    )

    QUALITY_CONTROL_AND_FILTERING(
        ch_raw_reads,
        ch_host_index,
        ch_host_index_prefix
    )

    MAG_CONSTRUCTION(
        QUALITY_CONTROL_AND_FILTERING.out.filtered_reads,
        ch_checkm2_db,
        ch_gunc_db
    )

    TAXONOMIC_CLASSIFICATION_AND_PHYLOGENOMICS(
        MAG_CONSTRUCTION.out.final_mags,
        ch_gtdbtk_db,
        ch_phylophlan_db,
        ch_phylophlan_config
    )

    GENE_PREDICTION_AND_FUNCTIONAL_ANNOTATION(
        MAG_CONSTRUCTION.out.final_mags,
        ch_genemark_home,
        ch_genemark_key,
        ch_eggnog_db,
        ch_interproscan_data
    )

    MAG_ABUNDANCE_ESTIMATION(
        MAG_CONSTRUCTION.out.final_mags,
        QUALITY_CONTROL_AND_FILTERING.out.filtered_reads
    )

    ch_global_reports = QUALITY_CONTROL_AND_FILTERING.out.raw_fastqc
        .mix(QUALITY_CONTROL_AND_FILTERING.out.fastp_json)
        .mix(QUALITY_CONTROL_AND_FILTERING.out.clean_fastqc)
        .mix(QUALITY_CONTROL_AND_FILTERING.out.bowtie2_logs)
        .mix(MAG_CONSTRUCTION.out.reports)
        .mix(MAG_CONSTRUCTION.out.logs)
        .mix(TAXONOMIC_CLASSIFICATION_AND_PHYLOGENOMICS.out.gtdbtk_summaries)
        .mix(TAXONOMIC_CLASSIFICATION_AND_PHYLOGENOMICS.out.logs)
        .mix(GENE_PREDICTION_AND_FUNCTIONAL_ANNOTATION.out.integration_summaries)
        .mix(GENE_PREDICTION_AND_FUNCTIONAL_ANNOTATION.out.logs)
        .mix(MAG_ABUNDANCE_ESTIMATION.out.abundance_wide)
        .mix(MAG_ABUNDANCE_ESTIMATION.out.abundance_long)
        .mix(MAG_ABUNDANCE_ESTIMATION.out.logs)

    GLOBAL_PROCESSING_EVALUATION(ch_global_reports, ch_multiqc_config)

    ch_versions = CHECK_SAMPLESHEET.out.versions
        .mix(QUALITY_CONTROL_AND_FILTERING.out.versions)
        .mix(MAG_CONSTRUCTION.out.versions)
        .mix(TAXONOMIC_CLASSIFICATION_AND_PHYLOGENOMICS.out.versions)
        .mix(GENE_PREDICTION_AND_FUNCTIONAL_ANNOTATION.out.versions)
        .mix(MAG_ABUNDANCE_ESTIMATION.out.versions)
        .mix(GLOBAL_PROCESSING_EVALUATION.out.versions)

    COLLECT_VERSIONS(ch_versions.collect(flat: false))

    emit:
    filtered_reads                  = QUALITY_CONTROL_AND_FILTERING.out.filtered_reads
    final_mags                      = MAG_CONSTRUCTION.out.final_mags
    final_catalog_provenance        = MAG_CONSTRUCTION.out.final_provenance
    final_catalog_quality           = MAG_CONSTRUCTION.out.final_quality
    final_checkm2                   = MAG_CONSTRUCTION.out.final_checkm2
    final_gunc                      = MAG_CONSTRUCTION.out.final_gunc
    final_species_representatives   = MAG_CONSTRUCTION.out.final_species_reps
    gtdbtk_summaries                = TAXONOMIC_CLASSIFICATION_AND_PHYLOGENOMICS.out.gtdbtk_summaries
    phylogenomic_tree               = TAXONOMIC_CLASSIFICATION_AND_PHYLOGENOMICS.out.phylophlan_tree
    functional_annotations         = GENE_PREDICTION_AND_FUNCTIONAL_ANNOTATION.out.integrated_annotations
    mag_abundance                   = MAG_ABUNDANCE_ESTIMATION.out.abundance_long
    multiqc_report                  = GLOBAL_PROCESSING_EVALUATION.out.report
    software_versions               = COLLECT_VERSIONS.out.table
}
