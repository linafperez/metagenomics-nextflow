process RESOLVE_SRA_PROJECT {
    tag "${project_accession}"
    label 'process_single'

    container 'python:3.12.11-slim-bookworm'
    conda "${moduleDir}/environment.yml"

    input:
    val project_accession
    val allowed_platforms
    val contact_email
    path resolver

    output:
    path 'sra_project_manifest.tsv', emit: run_manifest
    path 'sra_sample_manifest.tsv', emit: sample_manifest
    path 'sra_project_exclusions.tsv', emit: exclusions
    path 'sra_project_summary.json', emit: summary
    path 'sra_project_runinfo.csv', emit: runinfo
    tuple val("${task.process}"), val('python'), val('3.12.11'), emit: versions

    script:
    if (contact_email && !(contact_email.toString() ==~ /^[A-Za-z0-9._%+\-]+@[A-Za-z0-9.\-]+$/)) {
        error 'SRA contact email contains unsupported characters'
    }
    if (!(allowed_platforms.toString() ==~ /^[A-Za-z0-9_.\-]+(?:,[A-Za-z0-9_.\-]+)*$/)) {
        error 'SRA platform allowlist must contain only comma-separated platform names'
    }
    def email_arg = contact_email ? "--email '${contact_email}'" : ''

    """
    python3 "${resolver}" \
        "${project_accession}" \
        --output-dir . \
        --platforms "${allowed_platforms}" \
        ${email_arg} \
        --write-invalid-and-succeed
    """
}

process VALIDATE_SRA_PROJECT {
    tag "frozen project manifest"
    label 'process_single'

    container 'python:3.12.11-slim-bookworm'
    conda "${moduleDir}/environment.yml"

    input:
    path run_manifest
    path sample_manifest
    path exclusions
    path summary
    path runinfo
    path resolver

    output:
    path 'sra_project_manifest.validated', emit: sentinel
    tuple val("${task.process}"), val('sra_project_resolver'), val('1.0.0'), emit: versions

    script:
    """
    mkdir frozen_manifest
    cp "${run_manifest}" "${sample_manifest}" "${exclusions}" "${summary}" "${runinfo}" frozen_manifest/
    python3 "${resolver}" --validate-existing frozen_manifest
    printf 'validated\n' > sra_project_manifest.validated
    """
}
