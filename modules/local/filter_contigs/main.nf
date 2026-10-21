process FILTER_CONTIGS {
    tag "${meta.id}"
    label 'process_single'

    // Official Bash image includes BusyBox awk; Conda declares GNU awk.
    container 'bash:5.2.37'
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
    export LC_ALL=C
    if awk -v minimum_length=${minimum_length} \\
        -v stats_file="${prefix}.contig_filter.tsv" \\
        -f "${filter_script}" < "${contigs}" > "${prefix}.filtered.fa.partial"; then
        mv -- "${prefix}.filtered.fa.partial" "${prefix}.filtered.fa"
    else
        rm -f -- "${prefix}.filtered.fa.partial"
        exit 1
    fi
    """
    command

}
