process MEGAHIT {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/megahit:1.2.9--haf24da9_8'

    input:
    tuple val(meta), path(reads, arity: '2..*')

    output:
    tuple val(meta), path('*.contigs.fa'), emit: contigs
    tuple val(meta), path('*.megahit.log'), emit: log
    tuple val("${task.process}"), val('megahit'), val('1.2.9'), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id

    if (reads.size() % 2 != 0) {
        error "MEGAHIT requires an even number of paired-end read files"
    }

    def read_pairs    = reads.collate(2)
    def forward_reads = read_pairs.collect { pair -> pair[0] }.join(',')
    def reverse_reads = read_pairs.collect { pair -> pair[1] }.join(',')

    """
    echo "Starting MEGAHIT coassembly"

    megahit \
        -1 "${forward_reads}" \
        -2 "${reverse_reads}" \
        --presets meta-large \
        --min-contig-len 1000 \
        --num-cpu-threads ${task.cpus} \
        --mem-flag 1 \
        --memory 0.9 \
        ${args} \
        --out-dir "${prefix}.megahit" \
        > "${prefix}.megahit.log" 2>&1

    mv "${prefix}.megahit/final.contigs.fa" "${prefix}.contigs.fa"
    test -s "${prefix}.contigs.fa"
    rm -rf -- "${prefix}.megahit"
    """

}
