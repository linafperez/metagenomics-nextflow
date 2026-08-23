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
    def names = binFiles.collect { bin -> bin.simpleName }
    if (names.toSet().size() != names.size()) {
        error "CheckM2 input bin identifiers are not unique"
    }
    def links = binFiles.collect { bin ->
        "ln -s \"${bin}\" \"input_bins/${bin.simpleName}.fa\""
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

    stub:
    def prefix = task.ext.prefix ?: meta.id
    def binFiles = bins instanceof List ? bins : [bins]
    def rows = binFiles.withIndex().collect { bin, index ->
        "${bin.simpleName}\t${95.0 + (index % 3)}\t${1.0 + (index % 2)}\tNeural Network (Specific Model)\t11\t0.90\t2000\t300\t2400\t50.0\t8\tNone"
    }.join('\n')
    """
    mkdir -p "${prefix}.checkm2"
    printf 'Name\tCompleteness\tContamination\tCompleteness_Model_Used\tTranslation_Table_Used\tCoding_Density\tContig_N50\tAverage_Gene_Length\tGenome_Size\tGC_Content\tTotal_Coding_Sequences\tAdditional_Notes\n${rows}\n' > "${prefix}.checkm2/quality_report.tsv"
    cp "${prefix}.checkm2/quality_report.tsv" "${prefix}.checkm2.quality_report.tsv"
    printf 'CheckM2 stub evaluation completed\n' > "${prefix}.checkm2.log"
    """
}
