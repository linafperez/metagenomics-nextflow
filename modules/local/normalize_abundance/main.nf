process NORMALIZE_ABUNDANCE {
    tag "${meta.id}"
    label 'process_single'

    container 'python:3.12.11-slim-bookworm'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(abundance)
    path normalization_script

    output:
    tuple val(meta), path('*.mag_abundance.long.tsv'), emit: table
    tuple val("${task.process}"), val('normalize_coverm_abundance'), val('1.0.0'), emit: versions

    script:
    def prefix = task.ext.prefix ?: meta.id
    def command = """
    python3 "${normalization_script}" \\
        --input "${abundance}" \\
        --output "${prefix}.mag_abundance.long.tsv"
    """
    command

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    python3 "${normalization_script}" \\
        --input "${abundance}" \\
        --output "${prefix}.mag_abundance.long.tsv"
    """
}
