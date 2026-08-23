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

    stub:
    def prefix = task.ext.prefix ?: meta.id

    """
    set -euo pipefail

    mkdir -p "${prefix}.multiqc_data"

    printf '%s\n' \
        '<!doctype html>' \
        '<html lang="en">' \
        '<head><meta charset="utf-8"><title>Metagenomics MultiQC stub report</title></head>' \
        '<body><h1>Metagenomics MultiQC report</h1><p>Structured stub output.</p></body>' \
        '</html>' \
        > "${prefix}.multiqc.html"

    printf '%s\n' \
        '{' \
        '  "config_analysis_dir_abs": ["."],' \
        '  "config_report_title": "Shotgun metagenomics processing evaluation",' \
        '  "multiqc_version": "1.35",' \
        '  "report_saved_raw_data": {},' \
        '  "report_general_stats_data": [],' \
        '  "stub": true' \
        '}' \
        > "${prefix}.multiqc_data/multiqc_data.json"

    printf 'Sample\n' > "${prefix}.multiqc_data/multiqc_general_stats.txt"
    printf 'Module\tSource\nMultiQC\tstub\n' > "${prefix}.multiqc_data/multiqc_sources.txt"
    printf 'MultiQC\t1.35\n' > "${prefix}.multiqc_data/multiqc_versions.txt"
    printf 'MultiQC 1.35 structured stub report completed\n' > "${prefix}.multiqc.log"
    """
}
