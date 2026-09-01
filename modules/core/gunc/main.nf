process GUNC {
    tag "${meta.id}"
    label 'process_high'

    container 'quay.io/biocontainers/gunc:1.0.6--pyhdfd78af_1'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(bins, arity: '1..*')
    path database

    output:
    tuple val(meta), path('*.gunc'), optional: true, emit: results
    tuple val(meta), path('*.gunc.summary.tsv'), emit: summary
    tuple val(meta), path('*.gunc.log'), emit: log
    tuple val("${task.process}"), val('gunc'), val('1.0.6'), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id
    def keep_native_outputs = params.save_intermediates.toString().toBoolean()
    def binFiles = bins instanceof List ? bins : [bins]
    def names = binFiles.collect { bin -> bin.name.replaceFirst(/(?i)\.(fa|fna|fasta)(\.gz)?$/, '') }
    if (names.toSet().size() != names.size()) {
        error "GUNC input bin identifiers are not unique"
    }
    def links = binFiles.collect { bin ->
        def binId = bin.name.replaceFirst(/(?i)\.(fa|fna|fasta)(\.gz)?$/, '')
        "ln -s \"\$(readlink -f '${bin}')\" \"input_bins/${binId}.fa\""
    }.join('\n')

    """
    set -euo pipefail

    mkdir -p input_bins "${prefix}.gunc" "${prefix}.gunc.tmp"
    ${links}

    gunc run \\
        --input_dir input_bins \\
        --db_file "${database}" \\
        --out_dir "${prefix}.gunc" \\
        --temp_dir "${prefix}.gunc.tmp" \\
        --threads ${task.cpus} \\
        ${args} \\
        > "${prefix}.gunc.log" 2>&1

    summary=\$(find "${prefix}.gunc" -type f -name '*maxCSS_level.tsv' -print -quit)
    if [ -z "\$summary" ]; then
        echo 'GUNC did not create a maxCSS summary table' >&2
        exit 1
    fi
    cp "\$summary" "${prefix}.gunc.summary.tsv"
    test -s "${prefix}.gunc.summary.tsv"
    rm -rf -- input_bins "${prefix}.gunc.tmp"
    ${keep_native_outputs ? '' : "rm -rf -- '${prefix}.gunc'"}
    """

}
