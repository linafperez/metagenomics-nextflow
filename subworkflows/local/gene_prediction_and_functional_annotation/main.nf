#!/usr/bin/env nextflow

include { GENEMARKS2 } from '../../../modules/core/genemarks2/main'
include { EGGNOGMAPPER } from '../../../modules/core/eggnogmapper/main'
include { INTERPROSCAN } from '../../../modules/core/interproscan/main'
include { INTEGRATE_ANNOTATIONS } from '../../../modules/local/integrate_annotations/main'

workflow GENE_PREDICTION_AND_FUNCTIONAL_ANNOTATION {
    take:
    ch_final_mags
    ch_genemark_home
    ch_genemark_key
    ch_eggnog_db
    ch_interproscan_data

    main:
    ch_mags_per_genome = ch_final_mags.flatMap { meta, mags ->
        def mag_files = mags instanceof Collection ? mags as List : [mags]
        if (!mag_files) {
            error 'Functional annotation requires at least one final MAG'
        }

        def entries = mag_files.collect { mag ->
            def raw_id = mag.name.replaceFirst(/(?i)\.(fa|fna|fasta)(\.gz)?$/, '')
            def mag_id = raw_id.replaceAll(/[^A-Za-z0-9_.-]/, '_')
            if (!mag_id) {
                error "Could not derive a MAG identifier from '${mag.name}'"
            }
            [file: mag, id: mag_id]
        }

        if (entries.collect { entry -> entry.id }.toSet().size() != entries.size()) {
            error 'Final MAG identifiers are not unique after sanitization'
        }

        entries.collect { entry ->
            def mag_meta = meta + [
                id        : entry.id,
                mag_id    : entry.id,
                catalog_id: meta.catalog_id ?: meta.id
            ]
            tuple(mag_meta, entry.file)
        }
    }

    GENEMARKS2(
        ch_mags_per_genome,
        ch_genemark_home,
        ch_genemark_key
    )

    ch_predicted_proteins = GENEMARKS2.out.predictions.map { meta, proteins, _cds, _gff ->
        tuple(meta, proteins)
    }

    EGGNOGMAPPER(ch_predicted_proteins, ch_eggnog_db)
    INTERPROSCAN(ch_predicted_proteins, ch_interproscan_data)

    ch_gene_context = GENEMARKS2.out.predictions.map { meta, proteins, _cds, gff ->
        tuple(meta.mag_id, meta, proteins, gff)
    }

    ch_annotation_sources = EGGNOGMAPPER.out.annotations
        .map { meta, annotations -> tuple(meta.mag_id, annotations) }
        .join(
            INTERPROSCAN.out.tsv.map { meta, annotations ->
                tuple(meta.mag_id, annotations)
            }
        )

    ch_integration_input = ch_gene_context
        .join(ch_annotation_sources)
        .map { _mag_id, meta, proteins, gff, eggnog_annotations, interproscan_tsv ->
            tuple(meta, proteins, gff, eggnog_annotations, interproscan_tsv)
        }

    INTEGRATE_ANNOTATIONS(
        ch_integration_input,
        file("${projectDir}/bin/integrate_functional_annotations.py", checkIfExists: true)
    )

    ch_logs = GENEMARKS2.out.log
        .mix(EGGNOGMAPPER.out.log)
        .mix(INTERPROSCAN.out.log)
        .mix(INTEGRATE_ANNOTATIONS.out.log)

    ch_versions = GENEMARKS2.out.versions
        .mix(EGGNOGMAPPER.out.versions)
        .mix(INTERPROSCAN.out.versions)
        .mix(INTEGRATE_ANNOTATIONS.out.versions)

    emit:
    predictions             = GENEMARKS2.out.predictions
    eggnog_annotations      = EGGNOGMAPPER.out.annotations
    eggnog_seed_orthologs    = EGGNOGMAPPER.out.seed_orthologs
    eggnog_orthologs         = EGGNOGMAPPER.out.orthologs
    interproscan_tsv         = INTERPROSCAN.out.tsv
    interproscan_gff3        = INTERPROSCAN.out.gff3
    interproscan_json        = INTERPROSCAN.out.json
    integrated_annotations  = INTEGRATE_ANNOTATIONS.out.annotations
    integration_summaries   = INTEGRATE_ANNOTATIONS.out.summary
    logs                     = ch_logs
    versions                 = ch_versions
}
