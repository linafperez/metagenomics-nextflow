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
    sed "s|^data.directory=.*|data.directory=\${data_dir}|" \
        "\${properties}" \
        > "\${INTERPROSCAN_USER_DIR}/.interproscan-5/interproscan.properties"

    export _JAVA_OPTIONS="-Duser.home=\${INTERPROSCAN_USER_DIR}"
    export INTERPROSCAN_CONF="\${INTERPROSCAN_USER_DIR}/.interproscan-5/interproscan.properties"

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
    """

    stub:
    def prefix = task.ext.prefix ?: (meta.mag_id ?: meta.id)

    """
    set -euo pipefail

    protein_id=\$(awk '/^>/{sub(/^>/, ""); split(\$0, fields, /[[:space:]]+/); print fields[1]; exit}' "${proteins}")
    protein_id=\${protein_id:-stub_protein_1}

    printf '%s\tstub-md5\t20\tPfam\tPF00001\tStub protein domain\t1\t20\t1e-10\tT\t2026-01-01\tIPR000001\tStub InterPro entry\tGO:0003674\tKEGG:map00010\n' \
        "\${protein_id}" > "${prefix}.interproscan.tsv"

    printf '##gff-version 3\n%s\tInterProScan\tprotein_match\t1\t20\t1e-10\t+\t.\tID=match1;Name=PF00001;Dbxref=InterPro:IPR000001,GO:0003674\n' \
        "\${protein_id}" > "${prefix}.interproscan.gff3"
    printf '[{"xref":[{"id":"%s"}],"matches":[{"signature":{"accession":"PF00001","entry":{"accession":"IPR000001","description":"Stub InterPro entry"}}}]}]\n' \
        "\${protein_id}" > "${prefix}.interproscan.json"

    printf 'InterProScan stub annotation completed\n' > "${prefix}.interproscan.log"
    """
}
