#!/usr/bin/env nextflow

include { BINNING as MEGAHIT_BINNING_CORE } from '../binning/main'

workflow MEGAHIT_BINNING {
    take:
    ch_assembly
    ch_filtered_reads

    main:
    MEGAHIT_BINNING_CORE(ch_assembly, ch_filtered_reads)

    emit:
    filtered_contigs   = MEGAHIT_BINNING_CORE.out.filtered_contigs
    refined_bins       = MEGAHIT_BINNING_CORE.out.refined_bins
    comebin_bins       = MEGAHIT_BINNING_CORE.out.comebin_bins
    metabat2_bins      = MEGAHIT_BINNING_CORE.out.metabat2_bins
    semibin2_bins      = MEGAHIT_BINNING_CORE.out.semibin2_bins
    vamb_bins          = MEGAHIT_BINNING_CORE.out.vamb_bins
    coverm_bams        = MEGAHIT_BINNING_CORE.out.coverm_bams
    metabat_depth      = MEGAHIT_BINNING_CORE.out.metabat_depth
    vamb_abundance     = MEGAHIT_BINNING_CORE.out.vamb_abundance
    dastool_summary    = MEGAHIT_BINNING_CORE.out.dastool_summary
    dastool_evaluation = MEGAHIT_BINNING_CORE.out.dastool_evaluation
    reports            = MEGAHIT_BINNING_CORE.out.reports
    logs               = MEGAHIT_BINNING_CORE.out.logs
    versions           = MEGAHIT_BINNING_CORE.out.versions
}
