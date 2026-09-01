process FASTP {
    tag "${meta.id}"
    label 'process_medium'

    container 'community.wave.seqera.io/library/fastp:1.0.1--c8b87fe62dcc103c'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(reads, arity: '2')

    output:
    tuple val(meta), path('*_trimmed_*.fastq.gz', arity: '2'), emit: reads
    tuple val(meta), path('*.fastp.json'), emit: json
    tuple val(meta), path('*.fastp.html'), emit: html
    tuple val("${task.process}"), val('fastp'), val('1.0.1'), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id

    """
    fastp \\
        --in1 "${reads[0]}" \\
        --in2 "${reads[1]}" \\
        --out1 "${prefix}_trimmed_1.fastq.gz" \\
        --out2 "${prefix}_trimmed_2.fastq.gz" \\
        --json "${prefix}.fastp.json" \\
        --html "${prefix}.fastp.html" \\
        ${args}
    """

}
