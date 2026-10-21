#!/usr/bin/env nextflow

include { FETCH_SRA_RUNINFO } from '../../../modules/local/sra_runinfo/main'
include { SRA_ACQUISITION } from '../sra_acquisition/main'
include { QUALITY_CONTROL_AND_FILTERING } from '../quality_control_and_filtering/main'
include { METAGENOMICS_GLOBAL } from '../metagenomics_global/main'

workflow SRA_INPUT {
    main:
    def group_column = params.groupColumn?.toString()?.trim() ?: ''
    ch_selection = channel.fromPath(params.sraSamples, checkIfExists: true)
        .splitCsv(header: false, sep: '\t', strip: true)
        .toList()
        .map { rows -> [samples: InputMetadata.sraSamples(rows, group_column, params.sraProject, params.spadesCoassemblyMode == 'condition')] }

    ch_discovery_versions = channel.empty()
    if (params.sraRunInfo) {
        ch_runinfo = channel.value(file(params.sraRunInfo, checkIfExists: true))
    } else {
        FETCH_SRA_RUNINFO(
            ch_selection.map { _selection -> params.sraProject },
            channel.value(params.sraEmail ?: ''),
            channel.value(file("${projectDir}/bin/fetch_sra_runinfo.py", checkIfExists: true))
        )
        ch_runinfo = FETCH_SRA_RUNINFO.out.runinfo
        ch_discovery_versions = FETCH_SRA_RUNINFO.out.versions
    }

    ch_runinfo_rows = ch_runinfo.splitCsv(header: false, strip: true).toList().map { rows -> [rows: rows] }
    ch_resolution = ch_selection.combine(ch_runinfo_rows).map { selection, runinfo ->
        SraMetadata.resolve(selection.samples, runinfo.rows, params.sraProject, params.sraPlatforms)
    }
    ch_sample_metadata = ch_resolution.map { InputMetadata.metadataTsv(it.samples) }
        .collectFile(name: 'sample_metadata.tsv', storeDir: "${params.outdir}/pipeline_info", newLine: false)
    ch_exclusions = ch_resolution.map { SraMetadata.exclusionsTsv(it) }
        .collectFile(name: 'sra_exclusions.tsv', storeDir: "${params.outdir}/pipeline_info", newLine: false)

    // Load the whole cohort once; Nextflow schedules all selected runs and samples.
    ch_samples = ch_resolution.flatMap { it.samples }
    SRA_ACQUISITION(ch_samples, channel.value(file("${projectDir}/bin/validate_fastq_pair.awk", checkIfExists: true)))

    def host_index_prefix = file(params.host_bowtie2_index).name
    ch_host_index = channel.fromPath("${params.host_bowtie2_index}*.bt2*", checkIfExists: true)
        .collect().map { paths ->
            if (paths.size() != 6) error "Bowtie2 index '${params.host_bowtie2_index}' must contain exactly six files"
            paths.toList().sort { a, b -> a.name <=> b.name }
        }
    QUALITY_CONTROL_AND_FILTERING(SRA_ACQUISITION.out.reads, ch_host_index, channel.value(host_index_prefix))

    ch_reports = QUALITY_CONTROL_AND_FILTERING.out.raw_fastqc
        .mix(QUALITY_CONTROL_AND_FILTERING.out.fastp_json)
        .mix(QUALITY_CONTROL_AND_FILTERING.out.clean_fastqc)
        .mix(QUALITY_CONTROL_AND_FILTERING.out.bowtie2_logs)
    ch_versions = ch_discovery_versions.mix(SRA_ACQUISITION.out.versions)
        .mix(QUALITY_CONTROL_AND_FILTERING.out.versions)

    METAGENOMICS_GLOBAL(QUALITY_CONTROL_AND_FILTERING.out.filtered_reads, ch_reports, ch_versions, ch_sample_metadata)

    emit:
    filtered_reads = QUALITY_CONTROL_AND_FILTERING.out.filtered_reads
    sample_metadata = ch_sample_metadata
    exclusions = ch_exclusions
    final_mags = METAGENOMICS_GLOBAL.out.final_mags
    mag_abundance = METAGENOMICS_GLOBAL.out.mag_abundance
    multiqc_report = METAGENOMICS_GLOBAL.out.multiqc_report
    software_versions = METAGENOMICS_GLOBAL.out.software_versions
}
