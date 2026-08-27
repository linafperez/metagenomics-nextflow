process MULTIQC {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${task.ext.container ?: 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'}"

    input:
    tuple val(meta), path(reports, stageAs: '?/*', arity: '1..*')
    path multiqc_config

    output:
    tuple val(meta), path('*.multiqc.html'), emit: report
    tuple val(meta), path('*.multiqc_data'), emit: data
    tuple val(meta), path('*.multiqc.log'), emit: log
    tuple val("${task.process}"), val('multiqc'), val('1.35'), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id

    """
    set -euo pipefail

    echo "Building global MultiQC report"

    multiqc \
        . \
        --config "${multiqc_config}" \
        --outdir multiqc_output \
        --filename multiqc_report.html \
        --force \
        ${args} \
        2>&1 | tee "${prefix}.multiqc.log"

    mv multiqc_output/multiqc_report.html "${prefix}.multiqc.html"
    mv multiqc_output/multiqc_data "${prefix}.multiqc_data"
    """

}
