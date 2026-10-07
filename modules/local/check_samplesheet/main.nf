process CHECK_SAMPLESHEET {
    tag "${samplesheet.simpleName}"
    label 'process_single'

    container 'python:3.12.11-bookworm'
    conda "${moduleDir}/environment.yml"

    input:
    path samplesheet
    path validation_script
    path coassembly_helper
    val group_column

    output:
    path 'validated_samplesheet.csv', emit: csv
    path 'sample_metadata.tsv', emit: metadata
    tuple val("${task.process}"), val('python'), val('3.12'), emit: versions

    script:
    def validate_groups = params.spadesCoassemblyMode == 'condition' ? """
    python3 "${coassembly_helper}" validate-input --kind local \\
        --input "${samplesheet}" --group-column "${group_column}"
    """ : ''
    """
    ${validate_groups}
    python3 "${validation_script}" \\
        --input "${samplesheet}" \\
        --output validated_samplesheet.csv \\
        --metadata-output sample_metadata.tsv \\
        --group-column "${group_column}"
    """

}
