process CHECK_SAMPLESHEET {
    tag "${samplesheet.simpleName}"
    label 'process_single'

    container 'python:3.12.11-slim-bookworm'
    conda "${moduleDir}/environment.yml"

    input:
    path samplesheet
    path validation_script

    output:
    path 'validated_samplesheet.csv', emit: csv
    tuple val("${task.process}"), val('python'), val('3.12'), emit: versions

    script:
    """
    python3 "${validation_script}" \\
        --input "${samplesheet}" \\
        --output validated_samplesheet.csv
    """

}
