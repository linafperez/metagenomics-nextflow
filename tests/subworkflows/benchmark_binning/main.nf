#!/usr/bin/env nextflow

include { FILTER_CONTIGS } from '../../../modules/local/filter_contigs/main'
include { COVERM_CONTIG } from '../../../modules/core/coverm/main'
include { COMEBIN } from '../../../modules/core/comebin/main'
include { METABAT2 } from '../../../modules/core/metabat2/main'
include { SEMIBIN2 } from '../../../modules/core/semibin2/main'
include { VAMB } from '../../../modules/core/vamb/main'
include { DASTOOL } from '../../../modules/core/dastool/main'

workflow BENCHMARK_BINNING {
    take:
    ch_assembly
    ch_filtered_reads
    binner_strategy

    main:
    def strategy = binner_strategy.toString().toLowerCase()
    def valid_strategies = ['comebin', 'metabat2', 'semibin2', 'vamb', 'all']
    if (!valid_strategies.contains(strategy)) {
        error "Unsupported benchmark binner '${binner_strategy}'. Expected one of: ${valid_strategies.join(', ')}"
    }

    def pipeline_root = params.pipeline_root ?: projectDir
    ch_filter_script = channel.value(
        file("${pipeline_root}/bin/filter_fasta_by_length.py", checkIfExists: true)
    )

    FILTER_CONTIGS(
        ch_assembly,
        ch_filter_script,
        channel.value(params.min_contig_length)
    )

    ch_contigs = FILTER_CONTIGS.out.contigs
        .map { meta, contigs ->
            tuple(meta + [binning_strategy: strategy], contigs)
        }

    ch_read_collection = ch_filtered_reads
        .collect(flat: false)
        .map { records ->
            if (!records) {
                error 'Benchmark binning requires at least one paired-end sample'
            }
            def ordered = records.toList().sort { left, right ->
                left[0].id.toString() <=> right[0].id.toString()
            }
            ordered.each { record ->
                if (record.size() != 2 || !record[0].id || record[1].size() != 2) {
                    error 'Benchmark binning expects tuples of [meta, [read_1, read_2]]'
                }
            }
            tuple(
                ordered.collect { record -> record[0].id.toString() },
                ordered.collectMany { record -> record[1].toList() }
            )
        }

    ch_coverm_input = ch_contigs
        .combine(ch_read_collection)
        .map { meta, contigs, sample_ids, reads ->
            tuple(meta, contigs, sample_ids, reads)
        }

    COVERM_CONTIG(ch_coverm_input)

    ch_comebin_input = ch_contigs
        .map { meta, contigs -> tuple(meta.id, meta, contigs) }
        .join(COVERM_CONTIG.out.bams.map { meta, bams -> tuple(meta.id, bams) })
        .map { _id, meta, contigs, bams -> tuple(meta, contigs, bams) }

    ch_metabat2_input = ch_contigs
        .map { meta, contigs -> tuple(meta.id, meta, contigs) }
        .join(COVERM_CONTIG.out.metabat_depth.map { meta, depth -> tuple(meta.id, depth) })
        .map { _id, meta, contigs, depth -> tuple(meta, contigs, depth) }

    ch_semibin2_input = ch_contigs
        .map { meta, contigs -> tuple(meta.id, meta, contigs) }
        .join(COVERM_CONTIG.out.bams.map { meta, bams -> tuple(meta.id, bams) })
        .map { _id, meta, contigs, bams -> tuple(meta, contigs, bams) }

    ch_vamb_input = ch_contigs
        .map { meta, contigs -> tuple(meta.id, meta, contigs) }
        .join(COVERM_CONTIG.out.vamb_abundance.map { meta, abundance -> tuple(meta.id, abundance) })
        .map { _id, meta, contigs, abundance -> tuple(meta, contigs, abundance) }

    ch_reports = FILTER_CONTIGS.out.stats
        .mix(COVERM_CONTIG.out.metabat_depth)
        .mix(COVERM_CONTIG.out.vamb_abundance)
    ch_logs = COVERM_CONTIG.out.log
    ch_versions = FILTER_CONTIGS.out.versions.mix(COVERM_CONTIG.out.versions)
    ch_raw_bins = channel.empty()
    ch_binner_maps = channel.empty()
    ch_dastool_summary = channel.empty()

    if (strategy == 'comebin' || strategy == 'all') {
        COMEBIN(ch_comebin_input)
        ch_reports = ch_reports.mix(COMEBIN.out.contigs2bin)
        ch_logs = ch_logs.mix(COMEBIN.out.log)
        ch_versions = ch_versions.mix(COMEBIN.out.versions)
        ch_raw_bins = ch_raw_bins.mix(COMEBIN.out.bins)
        ch_binner_maps = ch_binner_maps.mix(COMEBIN.out.contigs2bin)
    }

    if (strategy == 'metabat2' || strategy == 'all') {
        METABAT2(ch_metabat2_input)
        ch_reports = ch_reports.mix(METABAT2.out.contigs2bin)
        ch_logs = ch_logs.mix(METABAT2.out.log)
        ch_versions = ch_versions.mix(METABAT2.out.versions)
        ch_raw_bins = ch_raw_bins.mix(METABAT2.out.bins)
        ch_binner_maps = ch_binner_maps.mix(METABAT2.out.contigs2bin)
    }

    if (strategy == 'semibin2' || strategy == 'all') {
        SEMIBIN2(ch_semibin2_input)
        ch_reports = ch_reports.mix(SEMIBIN2.out.contigs2bin)
        ch_logs = ch_logs.mix(SEMIBIN2.out.log)
        ch_versions = ch_versions.mix(SEMIBIN2.out.versions)
        ch_raw_bins = ch_raw_bins.mix(SEMIBIN2.out.bins)
        ch_binner_maps = ch_binner_maps.mix(SEMIBIN2.out.contigs2bin)
    }

    if (strategy == 'vamb' || strategy == 'all') {
        VAMB(ch_vamb_input)
        ch_reports = ch_reports.mix(VAMB.out.contigs2bin)
        ch_logs = ch_logs.mix(VAMB.out.log)
        ch_versions = ch_versions.mix(VAMB.out.versions)
        ch_raw_bins = ch_raw_bins.mix(VAMB.out.bins)
        ch_binner_maps = ch_binner_maps.mix(VAMB.out.contigs2bin)
    }

    if (strategy == 'all') {
        ch_dastool_input = ch_contigs
            .map { meta, contigs -> tuple(meta.id, meta, contigs) }
            .join(COMEBIN.out.contigs2bin.map { meta, mapping -> tuple(meta.id, mapping) })
            .join(METABAT2.out.contigs2bin.map { meta, mapping -> tuple(meta.id, mapping) })
            .join(SEMIBIN2.out.contigs2bin.map { meta, mapping -> tuple(meta.id, mapping) })
            .join(VAMB.out.contigs2bin.map { meta, mapping -> tuple(meta.id, mapping) })
            .map { _id, meta, contigs, comebin_map, metabat2_map, semibin2_map, vamb_map ->
                tuple(meta, contigs, comebin_map, metabat2_map, semibin2_map, vamb_map)
            }

        DASTOOL(ch_dastool_input)
        ch_refined_bins = DASTOOL.out.bins
        ch_dastool_summary = DASTOOL.out.summary
        ch_reports = ch_reports
            .mix(DASTOOL.out.summary)
            .mix(DASTOOL.out.evaluations)
        ch_logs = ch_logs.mix(DASTOOL.out.log)
        ch_versions = ch_versions.mix(DASTOOL.out.versions)
    }
    else if (strategy == 'comebin') {
        ch_refined_bins = COMEBIN.out.bins
    }
    else if (strategy == 'metabat2') {
        ch_refined_bins = METABAT2.out.bins
    }
    else if (strategy == 'semibin2') {
        ch_refined_bins = SEMIBIN2.out.bins
    }
    else {
        ch_refined_bins = VAMB.out.bins
    }

    emit:
    bins            = ch_refined_bins
    raw_bins        = ch_raw_bins
    binner_maps     = ch_binner_maps
    dastool_summary = ch_dastool_summary
    reports         = ch_reports
    logs            = ch_logs
    versions        = ch_versions
}
