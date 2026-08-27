process CHECKM2 {
    tag "${meta.id}"
    label 'process_high'

    container 'quay.io/biocontainers/checkm2:1.1.0--pyh7e72e81_1'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(bins, arity: '1..*')
    path database

    output:
    tuple val(meta), path('*.checkm2'), emit: results
    tuple val(meta), path('*.checkm2.quality_report.tsv'), emit: quality
    tuple val(meta), path('*.checkm2.log'), emit: log
    tuple val("${task.process}"), val('checkm2'), val('1.1.0'), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id
    def binFiles = bins instanceof List ? bins : [bins]
    def names = binFiles.collect { bin -> bin.name.replaceFirst(/(?i)\.(fa|fna|fasta)(\.gz)?$/, '') }
    if (names.toSet().size() != names.size()) {
        error "CheckM2 input bin identifiers are not unique"
    }
    def links = binFiles.collect { bin ->
        def binId = bin.name.replaceFirst(/(?i)\.(fa|fna|fasta)(\.gz)?$/, '')
        "ln -s \"\$(readlink -f '${bin}')\" \"input_bins/${binId}.fa\""
    }.join('\n')

    """
    mkdir -p input_bins
    ${links}

    database_path="${database}"
    if [ -d "\$database_path" ]; then
        database_path=\$(find -L "\$database_path" -type f -name '*.dmnd' -print -quit)
    fi
    if [ -z "\$database_path" ] || [ ! -f "\$database_path" ]; then
        echo 'CheckM2 DIAMOND database was not found' >&2
        exit 1
    fi

    checkm2 predict \\
        --input input_bins \\
        --output-directory "${prefix}.checkm2" \\
        --database_path "\$database_path" \\
        --threads ${task.cpus} \\
        -x fa \\
        ${args} \\
        > "${prefix}.checkm2.log" 2>&1

    cp "${prefix}.checkm2/quality_report.tsv" "${prefix}.checkm2.quality_report.tsv"
    """

}
