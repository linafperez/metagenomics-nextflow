process FILTER_CONTIGS {
    tag "${meta.id}"
    label 'process_single'

    container 'python:3.12.11-slim-bookworm'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(contigs)
    path filter_script
    val minimum_length

    output:
    tuple val(meta), path('*.filtered.fa'), emit: contigs
    tuple val(meta), path('*.contig_filter.tsv'), emit: stats
    tuple val("${task.process}"), val('filter_fasta_by_length'), val('1.0.0'), emit: versions

    script:
    def prefix = task.ext.prefix ?: meta.id
    def command = """
    python3 "${filter_script}" \\
        --input "${contigs}" \\
        --output "${prefix}.filtered.fa" \\
        --stats "${prefix}.contig_filter.tsv" \\
        --min-length ${minimum_length}
    """
    command

}
