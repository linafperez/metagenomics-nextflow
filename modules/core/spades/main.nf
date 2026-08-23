process SPADES {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/spades:4.2.0--h8d6e82b_2'

    input:
    tuple val(meta), path(reads, arity: '2..*')

    output:
    tuple val(meta), path('*.contigs.fa'), emit: contigs
    tuple val(meta), path('*.scaffolds.fa'), optional: true, emit: scaffolds
    tuple val(meta), path('*.assembly.gfa'), optional: true, emit: graph
    tuple val(meta), path('*.spades.log'), emit: log
    tuple val(meta), path('*.spades.params.txt'), emit: params
    tuple val("${task.process}"), val('spades'), val('4.2.0'), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id

    if (reads.size() % 2 != 0) {
        error "SPAdes requires an even number of paired-end read files"
    }

    def read_pairs  = reads.collate(2)
    def left_reads  = read_pairs.collect { pair -> '"' + pair[0] + '"' }.join(', ')
    def right_reads = read_pairs.collect { pair -> '"' + pair[1] + '"' }.join(', ')
    def dataset_yaml = """[
  {
    orientation: \"fr\",
    type: \"paired-end\",
    left reads: [${left_reads}],
    right reads: [${right_reads}]
  }
]"""
    def dataset_shell = dataset_yaml.replace("'", "'\"'\"'")
    def memory_gb     = task.memory ? Math.max(1, task.memory.toGiga() as int) : 1

    """
    echo "Starting metaSPAdes coassembly"

    printf '%s\\n' '${dataset_shell}' > "${prefix}.dataset.yaml"

    spades.py \
        --meta \
        --dataset "${prefix}.dataset.yaml" \
        --threads ${task.cpus} \
        --memory ${memory_gb} \
        ${args} \
        --output "${prefix}.spades" \
        > "${prefix}.spades.stdout.log" 2>&1

    mv "${prefix}.spades/contigs.fasta" "${prefix}.contigs.fa"
    mv "${prefix}.spades/spades.log" "${prefix}.spades.log"
    mv "${prefix}.spades/params.txt" "${prefix}.spades.params.txt"

    if [ -f "${prefix}.spades/scaffolds.fasta" ]; then
        mv "${prefix}.spades/scaffolds.fasta" "${prefix}.scaffolds.fa"
    fi

    if [ -f "${prefix}.spades/assembly_graph_with_scaffolds.gfa" ]; then
        mv "${prefix}.spades/assembly_graph_with_scaffolds.gfa" "${prefix}.assembly.gfa"
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id

    """
    {
        printf '>stub_spades_contig_1\n'
        count=0
        while [ "\$count" -lt 2400 ]; do
            printf 'TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA'
            count=\$((count + 64))
        done
        printf '\n'
    } > "${prefix}.contigs.fa"

    cp "${prefix}.contigs.fa" "${prefix}.scaffolds.fa"
    printf 'H\tVN:Z:1.0\nS\tstub_spades_contig_1\tTGCATGCA\n' > "${prefix}.assembly.gfa"
    printf 'SPAdes stub coassembly completed\n' > "${prefix}.spades.log"
    printf 'Command line: spades.py --meta --dataset stub.dataset.yaml\n' > "${prefix}.spades.params.txt"
    """
}
