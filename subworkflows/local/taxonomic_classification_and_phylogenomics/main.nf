#!/usr/bin/env nextflow

include { GTDBTK_CLASSIFY } from '../../../modules/core/gtdbtk/main'
include { PHYLOPHLAN } from '../../../modules/core/phylophlan/main'

workflow TAXONOMIC_CLASSIFICATION_AND_PHYLOGENOMICS {
    take:
    ch_final_mags
    ch_gtdbtk_db
    ch_phylophlan_db
    ch_phylophlan_config

    main:
    GTDBTK_CLASSIFY(ch_final_mags, ch_gtdbtk_db)

    PHYLOPHLAN(
        ch_final_mags,
        ch_phylophlan_db,
        ch_phylophlan_config
    )

    ch_gtdbtk_summaries = GTDBTK_CLASSIFY.out.bac120_summary
        .mix(GTDBTK_CLASSIFY.out.ar53_summary)

    ch_logs = GTDBTK_CLASSIFY.out.log
        .mix(PHYLOPHLAN.out.log)

    ch_versions = GTDBTK_CLASSIFY.out.versions
        .mix(PHYLOPHLAN.out.versions)
        .mix(PHYLOPHLAN.out.versions_iqtree)

    emit:
    gtdbtk_bac120_summary = GTDBTK_CLASSIFY.out.bac120_summary
    gtdbtk_ar53_summary   = GTDBTK_CLASSIFY.out.ar53_summary
    gtdbtk_summaries      = ch_gtdbtk_summaries
    gtdbtk_results        = GTDBTK_CLASSIFY.out.results
    phylophlan_tree       = PHYLOPHLAN.out.tree
    phylophlan_alignment  = PHYLOPHLAN.out.alignment
    phylophlan_results    = PHYLOPHLAN.out.results
    logs                  = ch_logs
    versions              = ch_versions
}
