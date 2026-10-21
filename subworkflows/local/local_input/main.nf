#!/usr/bin/env nextflow

include { QUALITY_CONTROL_AND_FILTERING } from '../quality_control_and_filtering/main'
include { METAGENOMICS_GLOBAL } from '../metagenomics_global/main'

workflow LOCAL_INPUT {
    main:
    def group_column = params.groupColumn?.toString()?.trim() ?: ''
    def samplesheet = file(params.input, checkIfExists: true)
    ch_cohort = channel.value(samplesheet).splitCsv(header: false, strip: true).toList()
        .map { rows ->
            InputMetadata.localSamples(rows, group_column, params.spadesCoassemblyMode == 'condition').collect { meta ->
                def reads = ['fastq_1', 'fastq_2'].collect { column ->
                    def source = meta.metadata[column]
                    if (!source) error "Sample ${meta.id} requires ${column}"
                    if (source.startsWith('~/')) source = System.getProperty('user.home') + source.substring(1)
                    def path = file(source.startsWith('/') ? source : samplesheet.parent.resolve(source).toString(), checkIfExists: true)
                    if (!(path.name.toLowerCase() ==~ /.+\.(fastq|fq)(\.gz)?$/) || !java.nio.file.Files.isRegularFile(path)) {
                        error "Invalid FASTQ file for ${meta.id}: ${path}"
                    }
                    path
                }
                if (reads[0] == reads[1]) error "Sample ${meta.id} uses the same file for both mates"
                tuple(meta, reads)
            }
        }
    ch_raw_reads = ch_cohort.flatMap { it }
    ch_sample_metadata = ch_cohort.map { records -> InputMetadata.metadataTsv(records.collect { it[0] }) }
        .collectFile(name: 'sample_metadata.tsv', storeDir: "${params.outdir}/pipeline_info", newLine: false)

    def host_index_prefix = file(params.host_bowtie2_index).name
    ch_host_index = channel.fromPath("${params.host_bowtie2_index}*.bt2*", checkIfExists: true)
        .collect().map { paths ->
            if (paths.size() != 6) error "Bowtie2 index '${params.host_bowtie2_index}' must contain exactly six files"
            paths.toList().sort { a, b -> a.name <=> b.name }
        }
    QUALITY_CONTROL_AND_FILTERING(ch_raw_reads, ch_host_index, channel.value(host_index_prefix))

    ch_preprocessing_reports = QUALITY_CONTROL_AND_FILTERING.out.raw_fastqc
        .mix(QUALITY_CONTROL_AND_FILTERING.out.fastp_json)
        .mix(QUALITY_CONTROL_AND_FILTERING.out.clean_fastqc)
        .mix(QUALITY_CONTROL_AND_FILTERING.out.bowtie2_logs)

    METAGENOMICS_GLOBAL(
        QUALITY_CONTROL_AND_FILTERING.out.filtered_reads,
        ch_preprocessing_reports,
        QUALITY_CONTROL_AND_FILTERING.out.versions,
        ch_sample_metadata
    )

    emit:
    filtered_reads = QUALITY_CONTROL_AND_FILTERING.out.filtered_reads
    sample_metadata = ch_sample_metadata
    final_mags = METAGENOMICS_GLOBAL.out.final_mags
    final_catalog_provenance = METAGENOMICS_GLOBAL.out.final_catalog_provenance
    final_catalog_quality = METAGENOMICS_GLOBAL.out.final_catalog_quality
    final_checkm2 = METAGENOMICS_GLOBAL.out.final_checkm2
    final_gunc = METAGENOMICS_GLOBAL.out.final_gunc
    final_species_representatives = METAGENOMICS_GLOBAL.out.final_species_representatives
    gtdbtk_summaries = METAGENOMICS_GLOBAL.out.gtdbtk_summaries
    phylogenomic_tree = METAGENOMICS_GLOBAL.out.phylogenomic_tree
    functional_annotations = METAGENOMICS_GLOBAL.out.functional_annotations
    mag_abundance = METAGENOMICS_GLOBAL.out.mag_abundance
    multiqc_report = METAGENOMICS_GLOBAL.out.multiqc_report
    software_versions = METAGENOMICS_GLOBAL.out.software_versions
}
