process EGGNOGMAPPER {
    tag "${meta.mag_id ?: meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${task.ext.container ?: 'quay.io/biocontainers/eggnog-mapper:2.1.13--pyhdfd78af_2'}"

    input:
    tuple val(meta), path(proteins)
    path eggnog_db

    output:
    tuple val(meta), path('*.emapper.annotations'), emit: annotations
    tuple val(meta), path('*.emapper.seed_orthologs'), emit: seed_orthologs
    tuple val(meta), path('*.emapper.orthologs'), emit: orthologs
    tuple val(meta), path('*.eggnogmapper.log'), emit: log
    tuple val("${task.process}"), val('eggnog-mapper'), val('2.1.13'), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: (meta.mag_id ?: meta.id)

    """
    set -euo pipefail

    database_dir=\$(readlink -f "${eggnog_db}")
    if [[ ! -s "\${database_dir}/eggnog.db" ]]; then
        echo "eggNOG 5.0.2 annotation database was not found" >&2
        exit 1
    fi
    if [[ ! -s "\${database_dir}/eggnog_proteins.dmnd" ]]; then
        echo "eggNOG 5.0.2 DIAMOND database was not found" >&2
        exit 1
    fi
    if [[ ! -s "\${database_dir}/eggnog.taxa.db" ]]; then
        echo "eggNOG 5.0.2 taxonomy database was not found" >&2
        exit 1
    fi

    echo "Running eggNOG-mapper functional annotation"
    mkdir -p eggnog_tmp
    export TMPDIR="\$PWD/eggnog_tmp"

    emapper.py \
        -i "${proteins}" \
        --itype proteins \
        -m diamond \
        --data_dir "\${database_dir}" \
        --output "${prefix}" \
        --output_dir . \
        --cpu ${task.cpus} \
        --report_orthologs \
        ${args} \
        > "${prefix}.eggnogmapper.log" 2>&1

    test -s "${prefix}.emapper.annotations"
    test -s "${prefix}.emapper.seed_orthologs"
    test -s "${prefix}.emapper.orthologs"
    rm -rf -- eggnog_tmp
    rm -f -- "${prefix}.emapper.hits"
    """

}
