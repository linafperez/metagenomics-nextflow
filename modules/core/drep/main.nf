process DREP {
    tag "${meta.id}:${stage}"
    label 'process_high'

    container 'quay.io/biocontainers/drep:3.6.2--pyhdfd78af_0'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(mags, arity: '1..*')
    path genome_info
    val ani
    val coverage
    val stage

    output:
    tuple val(meta), path('*.representatives/*.fa', arity: '1..*'), emit: representatives
    tuple val(meta), path('*.clusters.csv'), emit: clusters
    tuple val(meta), path('*.drep'), optional: true, emit: results
    tuple val(meta), path('*.drep.log'), emit: log
    tuple val("${task.process}"), val('drep'), val('3.6.2'), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def basePrefix = task.ext.prefix ?: meta.id
    def prefix = "${basePrefix}_${stage}"
    def keep_native_outputs = params.save_intermediates.toString().toBoolean()
    def magFiles = mags instanceof List ? mags : [mags]
    def links = magFiles.collect { mag ->
        def magId = mag.name.replaceFirst(/(?i)\.(fa|fna|fasta)(\.gz)?$/, '')
        "ln -s \"\$(readlink -f '${mag}')\" \"input_mags/${magId}.fa\""
    }.join('\n')

    """
    set -euo pipefail

    mkdir -p input_mags
    ${links}

    dRep dereplicate "${prefix}.drep" \\
        -g input_mags/*.fa \\
        --genomeInfo "${genome_info}" \\
        -p ${task.cpus} \\
        -comp 0 \\
        -con 100 \\
        -pa 0.90 \\
        -sa ${ani} \\
        -nc ${coverage} \\
        -cm larger \\
        --S_algorithm fastANI \\
        --multiround_primary_clustering \\
        --run_tertiary_clustering \\
        ${args} \\
        > "${prefix}.drep.log" 2>&1

    mkdir -p "${prefix}.representatives"
    cp "${prefix}.drep"/dereplicated_genomes/*.fa "${prefix}.representatives/"
    cp "${prefix}.drep/data_tables/Cdb.csv" "${prefix}.clusters.csv"
    test -n "\$(find "${prefix}.representatives" -type f -name '*.fa' -size +0c -print -quit)"
    test -s "${prefix}.clusters.csv"
    rm -rf -- input_mags
    ${keep_native_outputs ? '' : "rm -rf -- '${prefix}.drep'"}
    """

}
