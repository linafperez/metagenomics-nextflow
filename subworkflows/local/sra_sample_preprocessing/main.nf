#!/usr/bin/env nextflow

include { SRATOOLS_ACQUIRE } from '../../../modules/core/sratools/main'
include { PERSIST_SRA_CHECKPOINT } from '../../../modules/local/sra_preprocessing/main'
include { QUALITY_CONTROL_AND_FILTERING } from '../quality_control_and_filtering/main'

workflow SRA_SAMPLE_PREPROCESSING {
    main:
    def required = [
        sraManifest        : params.sraManifest,
        sraSampleId        : params.sraSampleId,
        sraCheckpointDir   : params.sraCheckpointDir,
        sraScratchDir      : params.sraScratchDir,
        sraCacheDir        : params.sraCacheDir,
        sraTempDir         : params.sraTempDir,
        host_bowtie2_index : params.host_bowtie2_index
    ]
    def missing = required.findAll { _name, value ->
        value == null || value.toString().trim().isEmpty()
    }.keySet().sort()
    if (missing) {
        error "SRA sample preprocessing is missing: ${missing.collect { name -> "--${name}" }.join(', ')}"
    }
    if (!(params.sraMaxSize.toString() ==~ /^(?:u|[1-9][0-9]*(?:[KMGT]B?)?)$/)) {
        error '--sraMaxSize must be u or a positive prefetch size'
    }

    ch_manifest = channel.value(file(params.sraManifest, checkIfExists: true))
    ch_sample_meta = channel
        .fromPath(params.sraManifest, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .filter { row -> row.sample_id == params.sraSampleId }
        .collect()
        .map { rows ->
            if (!rows) {
                error "SRA sample ${params.sraSampleId} is absent from the frozen manifest"
            }
            def groups = rows.collect { row -> row.group ?: '' }.toSet()
            def biosamples = rows.collect { row -> row.biosample_accession }.toSet()
            def hashes = rows.collect { row -> row.selection_file_sha256 }.toSet()
            if (groups.size() != 1 || biosamples.size() != 1 || hashes.size() != 1) {
                error "SRA sample ${params.sraSampleId} has contradictory frozen metadata"
            }
            [
                id: params.sraSampleId,
                single_end: false,
                biosample_accession: biosamples.first(),
                group: groups.first(),
                identity_source: 'BioSample',
                sample_order: rows.first().sample_order.toInteger(),
                run_accessions: rows.sort { left, right -> left.run_order.toInteger() <=> right.run_order.toInteger() }
                    .collect { row -> row.run_accession },
                selection_file_sha256: hashes.first()
            ]
        }

    ch_acquisition_helper = channel.value(
        file("${projectDir}/bin/acquire_sra_sample.py", checkIfExists: true)
    )
    ch_checkpoint_helper = channel.value(
        file("${projectDir}/bin/manage_sra_checkpoints.py", checkIfExists: true)
    )

    SRATOOLS_ACQUIRE(
        ch_sample_meta,
        ch_manifest,
        ch_acquisition_helper,
        channel.value(params.sraScratchDir),
        channel.value(params.sraCacheDir),
        channel.value(params.sraTempDir),
        channel.value(params.sraMaxSize)
    )

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

    QUALITY_CONTROL_AND_FILTERING(
        SRATOOLS_ACQUIRE.out.reads,
        ch_host_index,
        channel.value(host_index_prefix)
    )

    ch_reports = QUALITY_CONTROL_AND_FILTERING.out.raw_fastqc
        .map { _meta, paths -> paths }
        .mix(QUALITY_CONTROL_AND_FILTERING.out.fastp_json.map { _meta, path -> path })
        .mix(QUALITY_CONTROL_AND_FILTERING.out.clean_fastqc.map { _meta, paths -> paths })
        .mix(QUALITY_CONTROL_AND_FILTERING.out.bowtie2_logs.map { _meta, path -> path })
        .mix(SRATOOLS_ACQUIRE.out.versions)
        .flatten()
        .collect(flat: false)

    PERSIST_SRA_CHECKPOINT(
        QUALITY_CONTROL_AND_FILTERING.out.filtered_reads,
        ch_manifest,
        ch_reports,
        QUALITY_CONTROL_AND_FILTERING.out.versions.collect(flat: false),
        ch_checkpoint_helper,
        channel.value(params.sraCheckpointDir)
    )

    emit:
    record   = PERSIST_SRA_CHECKPOINT.out.record
    versions = PERSIST_SRA_CHECKPOINT.out.versions
}
