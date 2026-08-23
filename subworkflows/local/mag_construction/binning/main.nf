#!/usr/bin/env nextflow

include { FILTER_CONTIGS } from '../../../../modules/local/filter_contigs/main'
include { COVERM_CONTIG } from '../../../../modules/core/coverm/main'
include { COMEBIN } from '../../../../modules/core/comebin/main'
include { METABAT2 } from '../../../../modules/core/metabat2/main'
include { SEMIBIN2 } from '../../../../modules/core/semibin2/main'
include { VAMB } from '../../../../modules/core/vamb/main'
include { DASTOOL } from '../../../../modules/core/dastool/main'

workflow BINNING {
    take:
    ch_assembly
    ch_filtered_reads

    main:
    def pipeline_root = params.pipeline_root ?: projectDir
    ch_filter_script = channel.value(
        file("${pipeline_root}/bin/filter_fasta_by_length.py", checkIfExists: true)
    )

    FILTER_CONTIGS(
        ch_assembly,
        ch_filter_script,
        channel.value(params.min_contig_length)
    )

    ch_read_collection = ch_filtered_reads
        .collect(flat: false)
        .map { records ->
            if (!records) {
                error "Binning requires at least one paired-end sample"
            }
            def ordered = records.toList().sort { left, right ->
                left[0].id.toString() <=> right[0].id.toString()
            }
            ordered.each { record ->
                if (record.size() != 2 || !record[0].id || record[1].size() != 2) {
                    error "Binning expects tuples of [meta, [read_1, read_2]]"
                }
            }
            tuple(
                ordered.collect { record -> record[0].id.toString() },
                ordered.collectMany { record -> record[1].toList() }
            )
        }

    ch_coverm_input = FILTER_CONTIGS.out.contigs
        .combine(ch_read_collection)
        .map { meta, contigs, sample_ids, reads ->
            tuple(meta, contigs, sample_ids, reads)
        }

    COVERM_CONTIG(ch_coverm_input)

    ch_comebin_input = FILTER_CONTIGS.out.contigs
        .map { meta, contigs -> tuple(meta.id, meta, contigs) }
        .join(COVERM_CONTIG.out.bams.map { meta, bams -> tuple(meta.id, bams) })
        .map { _id, meta, contigs, bams -> tuple(meta, contigs, bams) }

    ch_metabat2_input = FILTER_CONTIGS.out.contigs
        .map { meta, contigs -> tuple(meta.id, meta, contigs) }
        .join(COVERM_CONTIG.out.metabat_depth.map { meta, depth -> tuple(meta.id, depth) })
        .map { _id, meta, contigs, depth -> tuple(meta, contigs, depth) }

    ch_semibin2_input = FILTER_CONTIGS.out.contigs
        .map { meta, contigs -> tuple(meta.id, meta, contigs) }
        .join(COVERM_CONTIG.out.bams.map { meta, bams -> tuple(meta.id, bams) })
        .map { _id, meta, contigs, bams -> tuple(meta, contigs, bams) }

    ch_vamb_input = FILTER_CONTIGS.out.contigs
        .map { meta, contigs -> tuple(meta.id, meta, contigs) }
        .join(COVERM_CONTIG.out.vamb_abundance.map { meta, abundance -> tuple(meta.id, abundance) })
        .map { _id, meta, contigs, abundance -> tuple(meta, contigs, abundance) }

    COMEBIN(ch_comebin_input)
    METABAT2(ch_metabat2_input)
    SEMIBIN2(ch_semibin2_input)
    VAMB(ch_vamb_input)

    ch_dastool_input = FILTER_CONTIGS.out.contigs
        .map { meta, contigs -> tuple(meta.id, meta, contigs) }
        .join(COMEBIN.out.contigs2bin.map { meta, mapping -> tuple(meta.id, mapping) })
        .join(METABAT2.out.contigs2bin.map { meta, mapping -> tuple(meta.id, mapping) })
        .join(SEMIBIN2.out.contigs2bin.map { meta, mapping -> tuple(meta.id, mapping) })
        .join(VAMB.out.contigs2bin.map { meta, mapping -> tuple(meta.id, mapping) })
        .map { _id, meta, contigs, comebin_map, metabat2_map, semibin2_map, vamb_map ->
            tuple(meta, contigs, comebin_map, metabat2_map, semibin2_map, vamb_map)
        }

    DASTOOL(ch_dastool_input)

    ch_reports = FILTER_CONTIGS.out.stats
        .mix(COVERM_CONTIG.out.metabat_depth)
        .mix(COVERM_CONTIG.out.vamb_abundance)
        .mix(COMEBIN.out.contigs2bin)
        .mix(METABAT2.out.contigs2bin)
        .mix(SEMIBIN2.out.contigs2bin)
        .mix(VAMB.out.contigs2bin)
        .mix(DASTOOL.out.summary)
        .mix(DASTOOL.out.evaluations)

    ch_logs = COVERM_CONTIG.out.log
        .mix(COMEBIN.out.log)
        .mix(METABAT2.out.log)
        .mix(SEMIBIN2.out.log)
        .mix(VAMB.out.log)
        .mix(DASTOOL.out.log)

    ch_versions = FILTER_CONTIGS.out.versions
        .mix(COVERM_CONTIG.out.versions)
        .mix(COMEBIN.out.versions)
        .mix(METABAT2.out.versions)
        .mix(SEMIBIN2.out.versions)
        .mix(VAMB.out.versions)
        .mix(DASTOOL.out.versions)

    emit:
    filtered_contigs  = FILTER_CONTIGS.out.contigs
    coverm_bams       = COVERM_CONTIG.out.bams
    metabat_depth     = COVERM_CONTIG.out.metabat_depth
    vamb_abundance    = COVERM_CONTIG.out.vamb_abundance
    comebin_bins      = COMEBIN.out.bins
    metabat2_bins     = METABAT2.out.bins
    semibin2_bins     = SEMIBIN2.out.bins
    vamb_bins         = VAMB.out.bins
    comebin_map       = COMEBIN.out.contigs2bin
    metabat2_map      = METABAT2.out.contigs2bin
    semibin2_map      = SEMIBIN2.out.contigs2bin
    vamb_map          = VAMB.out.contigs2bin
    refined_bins      = DASTOOL.out.bins
    dastool_summary   = DASTOOL.out.summary
    dastool_map       = DASTOOL.out.contigs2bin
    dastool_evaluation = DASTOOL.out.evaluations
    reports           = ch_reports
    logs              = ch_logs
    versions          = ch_versions
}
