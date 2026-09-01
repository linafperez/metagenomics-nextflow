#!/usr/bin/env nextflow

include { SRA_ACQUIRE; PERSIST_SRA_CHECKPOINT } from '../modules/local/sra_preprocessing/main'
include { QUALITY_CONTROL_AND_FILTERING } from '../subworkflows/local/quality_control_and_filtering/main'

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
    ch_acquisition_helper = channel.value(
        file("${projectDir}/bin/acquire_sra_sample.py", checkIfExists: true)
    )
    ch_checkpoint_helper = channel.value(
        file("${projectDir}/bin/manage_sra_checkpoints.py", checkIfExists: true)
    )

    SRA_ACQUIRE(
        channel.value(params.sraSampleId),
        ch_manifest,
        ch_acquisition_helper,
        channel.value(params.sraScratchDir),
        channel.value(params.sraCacheDir),
        channel.value(params.sraTempDir),
        channel.value(params.sraMaxSize)
    )

    ch_raw_reads = SRA_ACQUIRE.out.reads.map { sample_id, reads ->
        def ordered_reads = reads.toList().sort { left, right -> left.name <=> right.name }
        tuple([id: sample_id, single_end: false], ordered_reads)
    }

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
        ch_raw_reads,
        ch_host_index,
        channel.value(host_index_prefix)
    )

    ch_reports = QUALITY_CONTROL_AND_FILTERING.out.raw_fastqc
        .map { _meta, paths -> paths }
        .mix(QUALITY_CONTROL_AND_FILTERING.out.fastp_json.map { _meta, path -> path })
        .mix(QUALITY_CONTROL_AND_FILTERING.out.clean_fastqc.map { _meta, paths -> paths })
        .mix(QUALITY_CONTROL_AND_FILTERING.out.bowtie2_logs.map { _meta, path -> path })
        .mix(SRA_ACQUIRE.out.versions)
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
