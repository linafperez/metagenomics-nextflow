process MERGE_SRA_RUNS {
    tag "${meta.id}"
    label 'process_medium'

    container "${params.sraContainer}"
    containerOptions "${params.sraContainerOptions ?: ''}"
    conda "${moduleDir}/../../core/sratools/environment.yml"

    input:
    tuple val(meta), path(reads, arity: '2..*')

    output:
    tuple val(meta), path("${meta.id}_R*.fastq.gz", arity: 2), emit: reads
    tuple val("${task.process}"), val('pigz'), val('2.8'), emit: versions

    script:
    // Only runs of this BioSample are merged. Pair order was fixed in the channel.
    def commands = [0, 1].collect { mate ->
        def paths = (0..<meta.run_accessions.size()).collect { i -> "'${reads[2*i + mate]}'" }.join(' ')
        def target = "${meta.id}_R${mate + 1}.fastq.gz"
        meta.run_accessions.size() == 1
            ? "ln -s -- ${paths} '${target}'"
            : "pigz -dc ${paths} | pigz -p ${task.cpus} -c > '${target}'\npigz -t '${target}'"
    }.join('\n')
    """
    ${commands}
    """
}
