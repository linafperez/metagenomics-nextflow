process SELECT_HIGH_QUALITY_MAGS {
    tag "${meta.id}"
    label 'process_single'

    container 'python:3.12.11-slim-bookworm'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(bins, arity: '1..*'), path(quality_report)
    path selection_script
    val completeness_threshold
    val contamination_threshold

    output:
    tuple val(meta), path('hq_mags/*.fa', arity: '1..*'), emit: mags
    tuple val(meta), path('hq_mags.tsv'), emit: table
    tuple val(meta), path('genomeInfo.csv'), emit: genome_info
    tuple val("${task.process}"), val('select_high_quality_mags'), val('1.0.0'), emit: versions

    script:
    def binFiles = bins instanceof List ? bins : [bins]
    def binArgs = binFiles.collect { bin -> "\"${bin}\"" }.join(' ')
    def assembler = meta.assembler ?: meta.branch ?: 'unknown'
    def command = """
    python3 "${selection_script}" \\
        --bins ${binArgs} \\
        --quality "${quality_report}" \\
        --output-dir hq_mags \\
        --selected-table hq_mags.tsv \\
        --genome-info genomeInfo.csv \\
        --completeness ${completeness_threshold} \\
        --contamination ${contamination_threshold} \\
        --assembler "${assembler}"
    """
    command

}
