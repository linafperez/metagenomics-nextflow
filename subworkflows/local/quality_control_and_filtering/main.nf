#!/usr/bin/env nextflow

include { FASTQC as FASTQC_RAW } from '../../../modules/core/fastqc/main'
include { FASTP } from '../../../modules/core/fastp/main'
include { FASTQC as FASTQC_CLEAN } from '../../../modules/core/fastqc/main'
include { BOWTIE2 as BOWTIE2_HOST_REMOVAL } from '../../../modules/core/bowtie2/main'

workflow QUALITY_CONTROL_AND_FILTERING {
    take:
    ch_raw_reads
    ch_host_index
    ch_host_index_prefix

    main:
    FASTQC_RAW(ch_raw_reads)

    // Keying both channels by sample ID makes fastp wait for raw FastQC while
    // preserving the original read tuple. FastQC remains an assessment step.
    ch_raw_reads_keyed = ch_raw_reads.map { meta, reads ->
        tuple(meta.id, meta, reads)
    }

    ch_raw_fastqc_done = FASTQC_RAW.out.zip.map { meta, reports ->
        tuple(meta.id, reports)
    }

    ch_reads_after_raw_fastqc = ch_raw_reads_keyed
        .join(ch_raw_fastqc_done)
        .map { _sample_id, meta, reads, _reports ->
            tuple(meta, reads)
        }

    FASTP(ch_reads_after_raw_fastqc)
    FASTQC_CLEAN(FASTP.out.reads)

    // The second keyed join enforces the requested cleaned-read FastQC stage
    // before Bowtie2 without making the reusable FastQC module pass reads on.
    ch_clean_reads_keyed = FASTP.out.reads.map { meta, reads ->
        tuple(meta.id, meta, reads)
    }

    ch_clean_fastqc_done = FASTQC_CLEAN.out.zip.map { meta, reports ->
        tuple(meta.id, reports)
    }

    ch_reads_after_clean_fastqc = ch_clean_reads_keyed
        .join(ch_clean_fastqc_done)
        .map { _sample_id, meta, reads, _reports ->
            tuple(meta, reads)
        }

    BOWTIE2_HOST_REMOVAL(
        ch_reads_after_clean_fastqc,
        ch_host_index,
        ch_host_index_prefix
    )

    ch_versions = FASTQC_RAW.out.versions
        .mix(FASTP.out.versions)
        .mix(FASTQC_CLEAN.out.versions)
        .mix(BOWTIE2_HOST_REMOVAL.out.versions)

    emit:
    filtered_reads = BOWTIE2_HOST_REMOVAL.out.reads
    raw_fastqc     = FASTQC_RAW.out.zip
    fastp_json     = FASTP.out.json
    clean_fastqc   = FASTQC_CLEAN.out.zip
    bowtie2_logs   = BOWTIE2_HOST_REMOVAL.out.log
    versions       = ch_versions
}
