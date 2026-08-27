process FINALIZE_MAG_CATALOG {
    tag "${meta.id}"
    label 'process_single'

    container 'python:3.12.11-slim-bookworm'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(representatives, arity: '1..*'), path(provenance), path(quality)
    path finalize_script

    output:
    tuple val(meta), path('final_catalog/*.fa', arity: '1..*'), emit: mags
    tuple val(meta), path('final_catalog.provenance.tsv'), emit: provenance
    tuple val(meta), path('final_catalog.quality.tsv'), emit: quality
    tuple val("${task.process}"), val('finalize_mag_catalog'), val('1.0.0'), emit: versions

    script:
    def representativeFiles = representatives instanceof List ? representatives : [representatives]
    def representativeArgs = representativeFiles.collect { mag -> "\"${mag}\"" }.join(' ')
    def command = """
    python3 "${finalize_script}" \\
        --representatives ${representativeArgs} \\
        --provenance "${provenance}" \\
        --quality "${quality}" \\
        --output-dir final_catalog \\
        --output-provenance final_catalog.provenance.tsv \\
        --output-quality final_catalog.quality.tsv
    """
    command

}
