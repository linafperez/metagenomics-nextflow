process INTEGRATE_ANNOTATIONS {
    tag "${meta.mag_id ?: meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${task.ext.container ?: 'python:3.12.11-slim-bookworm'}"

    input:
    tuple val(meta), path(proteins), path(gff), path(eggnog_annotations), path(interproscan_tsv)
    path annotation_script

    output:
    tuple val(meta), path('*.functional_annotations.tsv'), emit: annotations
    tuple val(meta), path('*.functional_annotations.summary.json'), emit: summary
    tuple val(meta), path('*.annotation_integration.log'), emit: log
    tuple val("${task.process}"), val('integrate_functional_annotations'), val('1.0.0'), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: (meta.mag_id ?: meta.id)
    def mag_id = meta.mag_id ?: meta.id

    """
    set -euo pipefail

    python3 "${annotation_script}" \
        --proteins "${proteins}" \
        --gff "${gff}" \
        --eggnog "${eggnog_annotations}" \
        --interpro "${interproscan_tsv}" \
        --mag-id "${mag_id}" \
        --output "${prefix}.functional_annotations.tsv" \
        --summary "${prefix}.functional_annotations.summary.json" \
        ${args} \
        2>&1 | tee "${prefix}.annotation_integration.log"
    """

}
