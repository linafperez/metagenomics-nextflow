process GENEMARKS2 {
    tag "${meta.mag_id ?: meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${task.ext.container ?: 'docker.io/library/perl:5.44.0-slim-bookworm'}"

    input:
    tuple val(meta), path(mag)
    path genemark_home
    path genemark_key

    output:
    tuple val(meta), path('*.proteins.faa'), path('*.cds.fna'), path('*.genes.gff'), emit: predictions
    tuple val(meta), path('*.genemarks2.log'), emit: log
    tuple val("${task.process}"), val('genemarks2'), val('1.15'), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: (meta.mag_id ?: meta.id)
    def mag_id = meta.mag_id ?: meta.id

    """
    set -euo pipefail

    NF_MAG="${mag}"
    NF_CPUS="${task.cpus}"
    NF_PREFIX="${prefix}"
    NF_MAG_ID="${mag_id}"
    NF_GENEMARK_HOME="${genemark_home}"
    NF_GENEMARK_KEY="${genemark_key}"

    if [[ ! -d "\${NF_GENEMARK_HOME}" ]]; then
        echo "GeneMarkS-2 home must be a directory" >&2
        exit 1
    fi
    if [[ ! -s "\${NF_GENEMARK_KEY}" ]]; then
        echo "GeneMarkS-2 license key is missing or empty" >&2
        exit 1
    fi

    mkdir -p genemark_install genemark_user
    cp -aL "\${NF_GENEMARK_HOME}/." genemark_install/

    gms2_script=\$(find genemark_install -type f -name 'gms2.pl' -print -quit)
    if [[ -z "\${gms2_script}" ]]; then
        echo "gms2.pl was not found in the GeneMarkS-2 home" >&2
        exit 1
    fi
    gms2_script=\$(readlink -f "\${gms2_script}")
    gms2_dir=\$(dirname "\${gms2_script}")
    chmod u+x "\${gms2_script}" "\${gms2_dir}"/* 2>/dev/null || true

    install -m 600 "\${NF_GENEMARK_KEY}" genemark_user/.gmhmmp2_key
    cp genemark_user/.gmhmmp2_key genemark_user/.gm_key

    GENEMARK_USER_DIR="\$PWD/genemark_user"
    export GMHMMP2_KEY="\${GENEMARK_USER_DIR}/.gmhmmp2_key"
    export PATH="\${gms2_dir}:\$PATH"
    export LC_ALL=C
    export LANG=C

    gms2_usage=\$(perl "\${gms2_script}" 2>&1 || true)
    detected_version=\$(
        grep '^Version:' <<< "\${gms2_usage}" |
        head -n 1 |
        cut -d' ' -f2 || true
    )
    if [[ "\${detected_version}" != 1.15* ]]; then
        echo "GeneMarkS-2 1.15 is required; detected: \${detected_version:-unknown}" >&2
        exit 1
    fi

    mag_path=\$(readlink -f "\${NF_MAG}")
    output_root="\$PWD"

    echo "Running GeneMarkS-2 gene prediction"

    cd "\${GENEMARK_USER_DIR}"
    perl "\${gms2_script}" \
        --seq "\${mag_path}" \
        --genome-type bacteria \
        --threads "\${NF_CPUS}" \
        --format gff \
        --output "\${output_root}/\${NF_PREFIX}.raw.gff" \
        --fnn "\${output_root}/\${NF_PREFIX}.raw.fnn" \
        --faa "\${output_root}/\${NF_PREFIX}.raw.faa" \
        ${args} \
        > "\${output_root}/\${NF_PREFIX}.genemarks2.log" 2>&1
    cd "\${output_root}"

    awk -v MAG="\${NF_MAG_ID}" '
        /^>/ {
            gene=\$1
            sub(/^>/, "", gene)
            contig=(NF > 1 ? \$2 : "unknown_contig")
            details=\$0
            sub(/^>[^[:space:]]+[[:space:]]*/, "", details)
            printf ">%s|%s|%s", MAG, contig, gene
            if (details != "") printf " %s", details
            printf "\\n"
            next
        }
        { print }
    ' "\${NF_PREFIX}.raw.faa" > "\${NF_PREFIX}.proteins.faa"

    awk -v MAG="\${NF_MAG_ID}" '
        /^>/ {
            gene=\$1
            sub(/^>/, "", gene)
            contig=(NF > 1 ? \$2 : "unknown_contig")
            details=\$0
            sub(/^>[^[:space:]]+[[:space:]]*/, "", details)
            printf ">%s|%s|%s", MAG, contig, gene
            if (details != "") printf " %s", details
            printf "\\n"
            next
        }
        { print }
    ' "\${NF_PREFIX}.raw.fnn" > "\${NF_PREFIX}.cds.fna"

    awk -v MAG="\${NF_MAG_ID}" '
        BEGIN { FS=OFS="\\t"; feature_number=0 }
        /^#/ { print; next }
        NF < 9 { print; next }
        {
            feature_number++
            gene=""
            count=split(\$9, parts, ";")
            for (index=1; index<=count; index++) {
                item=parts[index]
                gsub(/^[[:space:]]+|[[:space:]]+\$/, "", item)
                if (item ~ /^gene_id[[:space:]=]/) {
                    sub(/^gene_id[[:space:]=]+/, "", item)
                    gsub(/^\"|\"\$/, "", item)
                    gene=item
                    break
                }
                if (item ~ /^ID=/) {
                    sub(/^ID=/, "", item)
                    gene=item
                }
            }
            if (gene == "") gene=feature_number
            protein_id=MAG "|" \$1 "|" gene
            if (\$9 !~ /;[[:space:]]*\$/) \$9=\$9 ";"
            \$9=\$9 "protein_id=" protein_id ";mag_id=" MAG
            print
        }
    ' "\${NF_PREFIX}.raw.gff" > "\${NF_PREFIX}.genes.gff"

    test -s "\${NF_PREFIX}.proteins.faa"
    test -s "\${NF_PREFIX}.cds.fna"
    test -s "\${NF_PREFIX}.genes.gff"
    rm -rf -- genemark_install genemark_user
    rm -f -- "\${NF_PREFIX}.raw.faa" "\${NF_PREFIX}.raw.fnn" "\${NF_PREFIX}.raw.gff"
    """

}
