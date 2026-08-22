process BOWTIE2 {
    tag "${meta.id}"
    label 'process_medium'

    container 'quay.io/biocontainers/bowtie2:2.5.4--he20e202_2'

    input:
    tuple val(meta), path(reads, arity: '2')
    path index, arity: '6'
    val index_prefix

    output:
    tuple val(meta), path('*_nonhost_*.fastq.gz', arity: '2'), emit: reads
    tuple val(meta), path('*.bowtie2.log'), emit: log
    tuple val("${task.process}"), val('bowtie2'),
        eval('bowtie2 --version | head -n 1 | sed "s/.*version //"'),
        emit: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id

    """
    bowtie2 \\
        -x "${index_prefix}" \\
        -1 "${reads[0]}" \\
        -2 "${reads[1]}" \\
        ${args} \\
        --threads ${task.cpus} \\
        --un-conc-gz "${prefix}_nonhost_%.fastq.gz" \\
        -S /dev/null \\
        2> "${prefix}.bowtie2.log"
    """
}
