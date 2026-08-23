#!/usr/bin/env nextflow

include { BINNING as SPADES_BINNING_CORE } from '../binning/main'

workflow SPADES_BINNING {
    take:
    ch_assembly
    ch_filtered_reads

    main:
    SPADES_BINNING_CORE(ch_assembly, ch_filtered_reads)

    emit:
    filtered_contigs   = SPADES_BINNING_CORE.out.filtered_contigs
    refined_bins       = SPADES_BINNING_CORE.out.refined_bins
    comebin_bins       = SPADES_BINNING_CORE.out.comebin_bins
    metabat2_bins      = SPADES_BINNING_CORE.out.metabat2_bins
    semibin2_bins      = SPADES_BINNING_CORE.out.semibin2_bins
    vamb_bins          = SPADES_BINNING_CORE.out.vamb_bins
    coverm_bams        = SPADES_BINNING_CORE.out.coverm_bams
    metabat_depth      = SPADES_BINNING_CORE.out.metabat_depth
    vamb_abundance     = SPADES_BINNING_CORE.out.vamb_abundance
    dastool_summary    = SPADES_BINNING_CORE.out.dastool_summary
    dastool_evaluation = SPADES_BINNING_CORE.out.dastool_evaluation
    reports            = SPADES_BINNING_CORE.out.reports
    logs               = SPADES_BINNING_CORE.out.logs
    versions           = SPADES_BINNING_CORE.out.versions
}
