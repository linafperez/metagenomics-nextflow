#!/usr/bin/env nextflow

include { SRATOOLS_ACQUIRE } from '../../../modules/core/sratools/main'
include { MERGE_SRA_RUNS } from '../../../modules/local/merge_sra_runs/main'

workflow SRA_ACQUISITION {
    take:
    ch_samples
    ch_pair_validator

    main:
    ch_runs = ch_samples.flatMap { meta -> meta.run_accessions.collect { run -> tuple(meta, run) } }
    SRATOOLS_ACQUIRE(ch_runs, ch_pair_validator)

    // groupKey emits each sample after its own expected run count is reached.
    ch_sample_runs = SRATOOLS_ACQUIRE.out.reads
        .map { meta, run, reads -> tuple(groupKey(meta, meta.run_accessions.size()), run, reads) }
        .groupTuple()
        .map { key, runs, pairs ->
            def meta = key.getGroupTarget()
            def ordered = meta.run_accessions.collectMany { run ->
                pairs[runs.indexOf(run)].toList().sort { a, b -> a.name <=> b.name }
            }
            tuple(meta, ordered)
        }
    MERGE_SRA_RUNS(ch_sample_runs)

    emit:
    reads = MERGE_SRA_RUNS.out.reads
    logs = SRATOOLS_ACQUIRE.out.log
    versions = SRATOOLS_ACQUIRE.out.versions.mix(MERGE_SRA_RUNS.out.versions)
}
