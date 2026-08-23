process FASTQC {
    tag "${meta.id}"
    label 'process_low'

    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(reads, arity: '2')

    output:
    tuple val(meta), path('*_fastqc.html', arity: '2'), emit: html
    tuple val(meta), path('*_fastqc.zip', arity: '2'), emit: zip
    tuple val("${task.process}"), val('fastqc'), val('0.12.1'), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id

    def renamed_reads = reads.withIndex().collect { read, index ->
        def name = read.name.toLowerCase()
        def extension

        if (name.endsWith('.fastq.gz')) {
            extension = '.fastq.gz'
        } else if (name.endsWith('.fq.gz')) {
            extension = '.fq.gz'
        } else if (name.endsWith('.fastq')) {
            extension = '.fastq'
        } else if (name.endsWith('.fq')) {
            extension = '.fq'
        } else {
            error "Unsupported FASTQ extension for '${read.name}'"
        }

        "${prefix}_${index + 1}${extension}"
    }

    def link_commands = reads.withIndex().collect { read, index ->
        "ln -s \"${read}\" \"${renamed_reads[index]}\""
    }.join('\n')

    """
    ${link_commands}

    fastqc \\
        --threads ${task.cpus} \\
        ${args} \\
        ${renamed_reads.join(' ')}
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    set -euo pipefail

    for mate in 1 2; do
        report="${prefix}_\${mate}_fastqc"
        mkdir -p "\${report}"
        printf '<html><body>FastQC stub</body></html>\n' > "\${report}/fastqc_report.html"
        printf '##FastQC\t0.12.1\n>>Basic Statistics\tpass\n#Measure\tValue\nFilename\t%s.fastq.gz\nFile type\tConventional base calls\nEncoding\tSanger / Illumina 1.9\nTotal Sequences\t1\nSequences flagged as poor quality\t0\nSequence length\t50\n%%GC\t50\n>>END_MODULE\n' \
            "${prefix}_\${mate}" > "\${report}/fastqc_data.txt"
        printf 'PASS\tBasic Statistics\t%s.fastq.gz\n' "${prefix}_\${mate}" > "\${report}/summary.txt"
        cp "\${report}/fastqc_report.html" "\${report}.html"
        python3 -m zipfile -c "\${report}.zip" "\${report}"
    done
    """
}
