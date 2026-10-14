process INTERPROSCAN {
    tag "${meta.mag_id ?: meta.id}"
    label 'process_high_memory'

    conda "${moduleDir}/environment.yml"
    container "${task.ext.container ?: 'quay.io/biocontainers/interproscan:5.59_91.0--hec16e2b_1'}"

    input:
    tuple val(meta), path(proteins)
    path interproscan_data

    output:
    tuple val(meta), path('*.interproscan.tsv'), emit: tsv
    tuple val(meta), path('*.interproscan.gff3'), emit: gff3
    tuple val(meta), path('*.interproscan.json'), emit: json
    tuple val(meta), path('*.interproscan.log'), emit: log
    tuple val("${task.process}"), val('interproscan'), val('5.59-91.0'), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: (meta.mag_id ?: meta.id)
    def fingerprintscan_path = task.ext.fingerprintscan_path ?: ''
    def cath_resolve_hits_path = task.ext.cath_resolve_hits_path ?: ''

    """
    set -euo pipefail

    data_dir=\$(readlink -f "${interproscan_data}")
    if [[ -d "\${data_dir}/data" ]]; then
        data_dir="\${data_dir}/data"
    fi
    if [[ ! -s "\${data_dir}/pfam/35.0/pfam_a.hmm" ]]; then
        echo "InterProScan 5.59-91.0 data were not found" >&2
        exit 1
    fi

    interproscan_script=\$(command -v interproscan.sh)
    interproscan_dir=\$(dirname "\$(readlink -f "\${interproscan_script}")")
    properties="\${interproscan_dir}/interproscan.properties"
    if [[ ! -s "\${properties}" ]]; then
        echo "InterProScan properties file was not found" >&2
        exit 1
    fi

    INTERPROSCAN_USER_DIR="\$PWD/interproscan_user"
    mkdir -p "\${INTERPROSCAN_USER_DIR}/.interproscan-5" tmp

    custom_properties="\${INTERPROSCAN_USER_DIR}/.interproscan-5/interproscan.properties"

    sed "s|^data.directory=.*|data.directory=\${data_dir}|" \
        "\${properties}" \
        > "\${custom_properties}"

    fingerprintscan_path="${fingerprintscan_path}"

    if [[ -n "\${fingerprintscan_path}" ]]; then
        fingerprintscan_path=\$(readlink -f "\${fingerprintscan_path}")

        if [[ ! -x "\${fingerprintscan_path}" ]]; then
            echo "Configured PRINTS fingerPRINTScan executable is not executable: \${fingerprintscan_path}" >&2
            exit 1
        fi

        if grep -q '^binary.fingerprintscan.path=' "\${custom_properties}"; then
            sed -i \
                "s|^binary.fingerprintscan.path=.*|binary.fingerprintscan.path=\${fingerprintscan_path}|" \
                "\${custom_properties}"
        else
            printf 'binary.fingerprintscan.path=%s\n' \
                "\${fingerprintscan_path}" >> "\${custom_properties}"
        fi
    fi

    cath_resolve_hits_path="${cath_resolve_hits_path}"

    if [[ -n "\${cath_resolve_hits_path}" ]]; then
        cath_resolve_hits_path=\$(readlink -f "\${cath_resolve_hits_path}")
        if [[ ! -x "\${cath_resolve_hits_path}" ]]; then
            echo "Configured CATH cath-resolve-hits executable is not executable: \${cath_resolve_hits_path}" >&2
            exit 1
        fi

        if grep -q '^cath.resolve.hits.path=' "\${custom_properties}"; then
            sed -i \
                "s|^cath.resolve.hits.path=.*|cath.resolve.hits.path=\${cath_resolve_hits_path}|" \
                "\${custom_properties}"
        else
            printf 'cath.resolve.hits.path=%s\n' \
                "\${cath_resolve_hits_path}" >> "\${custom_properties}"
        fi
    fi

    export _JAVA_OPTIONS="-Duser.home=\${INTERPROSCAN_USER_DIR}"
    export INTERPROSCAN_CONF="\${custom_properties}"

    echo "Running InterProScan functional annotation"

    interproscan.sh \
        -i "${proteins}" \
        -t p \
        -f TSV,GFF3,JSON \
        -b "${prefix}.interproscan" \
        -cpu ${task.cpus} \
        -T "\$PWD/tmp" \
        -dp \
        -goterms \
        -pa \
        ${args} \
        > "${prefix}.interproscan.log" 2>&1

    test -f "${prefix}.interproscan.tsv"
    test -s "${prefix}.interproscan.gff3"
    test -s "${prefix}.interproscan.json"
    rm -rf -- tmp interproscan_user
    """

}
