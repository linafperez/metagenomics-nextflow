process BOWTIE2 {
    tag "${meta.id}"
    label 'process_medium'

    container 'quay.io/biocontainers/bowtie2:2.5.4--he20e202_2'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(reads, arity: '2')
    path index, arity: '6'
    val index_prefix

    output:
    tuple val(meta), path('*_nonhost_*.fastq.gz', arity: '2'), emit: reads
    tuple val(meta), path('*.bowtie2.log'), emit: log
    tuple val("${task.process}"), val('bowtie2'), val('2.5.4'), emit: versions

    when:
    task.ext.when == null || task.ext.when

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

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    printf '@stub/1\nACGTACGTACGT\n+\nFFFFFFFFFFFF\n' | gzip -c > "${prefix}_nonhost_1.fastq.gz"
    printf '@stub/2\nTGCATGCATGCA\n+\nFFFFFFFFFFFF\n' | gzip -c > "${prefix}_nonhost_2.fastq.gz"
    printf '2 reads; of these:\n  2 (100.00%%) were paired; of these:\n    2 (100.00%%) aligned concordantly 0 times\n0.00%% overall alignment rate\n' > "${prefix}.bowtie2.log"
    """
}
