process GUNC {
    tag "${meta.id}"
    label 'process_high'

    container 'quay.io/biocontainers/gunc:1.0.6--pyhdfd78af_1'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(bins, arity: '1..*')
    path database

    output:
    tuple val(meta), path('*.gunc'), emit: results
    tuple val(meta), path('*.gunc.summary.tsv'), emit: summary
    tuple val(meta), path('*.gunc.log'), emit: log
    tuple val("${task.process}"), val('gunc'), val('1.0.6'), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id
    def binFiles = bins instanceof List ? bins : [bins]
    def names = binFiles.collect { bin -> bin.simpleName }
    if (names.toSet().size() != names.size()) {
        error "GUNC input bin identifiers are not unique"
    }
    def links = binFiles.collect { bin ->
        "ln -s \"${bin}\" \"input_bins/${bin.simpleName}.fa\""
    }.join('\n')

    """
    mkdir -p input_bins
    ${links}

    gunc run \\
        --input_dir input_bins \\
        --db_file "${database}" \\
        --out_dir "${prefix}.gunc" \\
        --threads ${task.cpus} \\
        ${args} \\
        > "${prefix}.gunc.log" 2>&1

    summary=\$(find "${prefix}.gunc" -type f -name '*maxCSS_level.tsv' -print -quit)
    if [ -z "\$summary" ]; then
        echo 'GUNC did not create a maxCSS summary table' >&2
        exit 1
    fi
    cp "\$summary" "${prefix}.gunc.summary.tsv"
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    def binFiles = bins instanceof List ? bins : [bins]
    def rows = binFiles.withIndex().collect { bin, index ->
        "${bin.simpleName}\t${100 + index}\t0.02\t0.01\t0.05\tgenus\t0.90\tTrue"
    }.join('\n')
    """
    mkdir -p "${prefix}.gunc"
    printf 'genome\tn_genes_called\tcontamination_portion\tn_effective_surplus_clades\tclade_separation_score\ttaxonomic_level\tproportion_genes_retained_in_major_clades\tpass.GUNC\n${rows}\n' > "${prefix}.gunc/GUNC.GTDB.maxCSS_level.tsv"
    cp "${prefix}.gunc/GUNC.GTDB.maxCSS_level.tsv" "${prefix}.gunc.summary.tsv"
    printf 'GUNC stub evaluation completed\n' > "${prefix}.gunc.log"
    """
}
