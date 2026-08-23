#!/usr/bin/env nextflow

include { COVERM_GENOME as COVERM_MAG_ABUNDANCE } from '../../../modules/core/coverm/main'
include { NORMALIZE_ABUNDANCE } from '../../../modules/local/normalize_abundance/main'

workflow MAG_ABUNDANCE_ESTIMATION {
    take:
    ch_final_mags
    ch_filtered_reads

    main:
    def pipeline_root = params.pipeline_root ?: projectDir

    ch_read_collection = ch_filtered_reads
        .collect(flat: false)
        .map { records ->
            if (!records) {
                error "MAG abundance estimation requires at least one paired-end sample"
            }
            def ordered = records.toList().sort { left, right ->
                left[0].id.toString() <=> right[0].id.toString()
            }
            tuple(
                ordered.collect { record -> record[0].id.toString() },
                ordered.collectMany { record -> record[1].toList() }
            )
        }

    ch_coverm_input = ch_final_mags
        .combine(ch_read_collection)
        .map { meta, mags, sample_ids, reads -> tuple(meta, mags, sample_ids, reads) }

    COVERM_MAG_ABUNDANCE(ch_coverm_input)

    ch_normalization_script = channel.value(
        file("${pipeline_root}/bin/normalize_coverm_abundance.py", checkIfExists: true)
    )
    NORMALIZE_ABUNDANCE(COVERM_MAG_ABUNDANCE.out.abundance, ch_normalization_script)

    ch_versions = COVERM_MAG_ABUNDANCE.out.versions
        .mix(NORMALIZE_ABUNDANCE.out.versions)

    emit:
    abundance_wide = COVERM_MAG_ABUNDANCE.out.abundance
    abundance_long = NORMALIZE_ABUNDANCE.out.table
    logs           = COVERM_MAG_ABUNDANCE.out.log
    versions       = ch_versions
}
