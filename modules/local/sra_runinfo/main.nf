process FETCH_SRA_RUNINFO {
    tag "${project_accession}"
    label 'process_single'

    container 'python:3.12.11-bookworm'
    conda "${moduleDir}/environment.yml"

    input:
    val project_accession
    val contact_email
    path fetcher

    output:
    path 'sra_runinfo.csv', emit: runinfo
    tuple val("${task.process}"), val('python'), val('3.12.11'), emit: versions

    script:
    def email_arg = contact_email ? "--email '${contact_email}'" : ''
    """
    python3 "${fetcher}" '${project_accession}' --output sra_runinfo.csv ${email_arg}
    """
}
