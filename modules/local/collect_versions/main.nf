process COLLECT_VERSIONS {
    tag 'software versions'
    label 'process_single'

    container 'python:3.12.11-slim-bookworm'
    conda "${moduleDir}/environment.yml"

    input:
    val version_records

    output:
    path 'software_versions.tsv', emit: table
    tuple val("${task.process}"), val('collect_versions'), val('1.0.0'), emit: versions

    script:
    def unique_records = (version_records ?: [])
        .collect { record -> record instanceof Collection ? record as List : [record] }
        .findAll { record -> record.size() >= 3 }
        .collect { record -> [record[0].toString(), record[1].toString(), record[2].toString()] }
        .unique { record -> record.join('\t') }
        .sort { left, right ->
            def tool_order = left[1] <=> right[1]
            tool_order != 0 ? tool_order : left[0] <=> right[0]
        }
    def rows = unique_records.collect { record -> record.join('\t') }.join('\n')
    def body = rows ? "${rows}\n" : ''

    """
    printf 'process\ttool\tversion\n%b' '${body}' > software_versions.tsv
    """

    stub:
    def unique_records = (version_records ?: [])
        .collect { record -> record instanceof Collection ? record as List : [record] }
        .findAll { record -> record.size() >= 3 }
        .collect { record -> [record[0].toString(), record[1].toString(), record[2].toString()] }
        .unique { record -> record.join('\t') }
        .sort { left, right -> left[1] <=> right[1] }
    def rows = unique_records.collect { record -> record.join('\t') }.join('\n')
    def body = rows ? "${rows}\n" : ''

    """
    printf 'process\ttool\tversion\n%b' '${body}' > software_versions.tsv
    """
}
