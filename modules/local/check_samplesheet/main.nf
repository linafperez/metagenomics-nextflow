process CHECK_SAMPLESHEET {
    tag "${samplesheet.simpleName}"
    label 'process_single'

    input:
    path samplesheet
    path validation_script

    output:
    path 'validated_samplesheet.csv', emit: csv
    tuple val("${task.process}"), val('python'),
        eval('python3 --version 2>&1 | awk \'{print $2}\''),
        emit: versions

    script:
    """
    python3 "${validation_script}" \\
        --input "${samplesheet}" \\
        --output validated_samplesheet.csv
    """
}
