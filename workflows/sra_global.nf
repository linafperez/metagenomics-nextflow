#!/usr/bin/env nextflow

include { METAGENOMICS_GLOBAL } from '../subworkflows/local/metagenomics_global/main'
include { FINALIZE_SRA_GLOBAL_RUN } from '../modules/local/sra_preprocessing/main'

workflow SRA_GLOBAL {
    main:
    if (!params.sraCheckpointManifest || !params.sraProject) {
        error 'SRA global analysis requires --sraCheckpointManifest and --sraProject'
    }

    ch_checkpoint_rows = channel
        .fromPath(params.sraCheckpointManifest, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')

    ch_filtered_reads = ch_checkpoint_rows.map { row ->
        def meta = [
            id: row.sample_id,
            single_end: false,
            biosample_accession: row.biosample_accession,
            identity_source: row.identity_source,
            sample_order: row.sample_order.toInteger(),
            run_accessions: row.run_accessions.tokenize(';')
        ]
        tuple(meta, [
            file(row.read_1, checkIfExists: true),
            file(row.read_2, checkIfExists: true)
        ])
    }

    ch_persisted_reports = ch_checkpoint_rows
        .flatMap { row ->
            new groovy.json.JsonSlurper().parseText(row.reports_json).findAll { path ->
                !path.toString().endsWith('_versions.tsv')
            }
        }
        .map { report -> file(report, checkIfExists: true) }

    ch_preprocessing_versions = ch_checkpoint_rows
        .flatMap { row ->
            new groovy.json.JsonSlurper().parseText(row.reports_json).findAll { path ->
                path.toString().endsWith('_versions.tsv') || path.toString().endsWith('preprocessing_versions.tsv')
            }
        }
        .map { version_file -> file(version_file, checkIfExists: true) }
        .splitCsv(header: false, sep: '\t')
        .map { fields -> tuple(fields[0], fields[1], fields[2]) }

    ch_static_versions = channel.of(
        tuple('RESOLVE_SRA_PROJECT', 'python', '3.12.11'),
        tuple('VALIDATE_SRA_PROJECT', 'sra_project_resolver', '1.0.0'),
        tuple('CHECK_SRA_CHECKPOINTS', 'sra_checkpoint_manager', '1.0.0')
    )

    METAGENOMICS_GLOBAL(
        ch_filtered_reads,
        ch_persisted_reports,
        ch_preprocessing_versions.mix(ch_static_versions)
    )

    ch_multiqc_html = METAGENOMICS_GLOBAL.out.multiqc_report.map { _meta, report ->
        report
    }
    ch_abundance_long = METAGENOMICS_GLOBAL.out.mag_abundance.map { _meta, table ->
        table
    }
    def durable_multiqc_path = file(
        "${params.outdir}/06_global_processing_evaluation/global_processing_evaluation.multiqc.html"
    ).toString()
    def durable_versions_path = file("${params.outdir}/pipeline_info/software_versions.tsv").toString()
    def durable_abundance_path = file(
        "${params.outdir}/05_mag_abundance_estimation/final_catalog.mag_abundance.long.tsv"
    ).toString()
    def durable_checkpoint_manifest_path = file(
        params.sraCheckpointManifest,
        checkIfExists: true
    ).toString()

    FINALIZE_SRA_GLOBAL_RUN(
        channel.value(params.sraProject),
        channel.value(file(params.sraCheckpointManifest, checkIfExists: true)),
        ch_multiqc_html,
        METAGENOMICS_GLOBAL.out.software_versions,
        ch_abundance_long,
        channel.value(durable_checkpoint_manifest_path),
        channel.value(durable_multiqc_path),
        channel.value(durable_versions_path),
        channel.value(durable_abundance_path)
    )

    emit:
    final_mags                    = METAGENOMICS_GLOBAL.out.final_mags
    final_catalog_provenance      = METAGENOMICS_GLOBAL.out.final_catalog_provenance
    final_catalog_quality         = METAGENOMICS_GLOBAL.out.final_catalog_quality
    final_checkm2                 = METAGENOMICS_GLOBAL.out.final_checkm2
    final_gunc                    = METAGENOMICS_GLOBAL.out.final_gunc
    final_species_representatives = METAGENOMICS_GLOBAL.out.final_species_representatives
    gtdbtk_summaries              = METAGENOMICS_GLOBAL.out.gtdbtk_summaries
    phylogenomic_tree             = METAGENOMICS_GLOBAL.out.phylogenomic_tree
    functional_annotations       = METAGENOMICS_GLOBAL.out.functional_annotations
    mag_abundance                 = METAGENOMICS_GLOBAL.out.mag_abundance
    multiqc_report                = METAGENOMICS_GLOBAL.out.multiqc_report
    software_versions             = METAGENOMICS_GLOBAL.out.software_versions
    success_marker                = FINALIZE_SRA_GLOBAL_RUN.out.marker
}
