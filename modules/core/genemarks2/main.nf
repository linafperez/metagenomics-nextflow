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

    if [[ ! -d "${genemark_home}" ]]; then
        echo "GeneMarkS-2 home must be a directory" >&2
        exit 1
    fi
    if [[ ! -s "${genemark_key}" ]]; then
        echo "GeneMarkS-2 license key is missing or empty" >&2
        exit 1
    fi

    mkdir -p genemark_install genemark_user
    cp -aL "${genemark_home}/." genemark_install/

    gms2_script=\$(find genemark_install -type f -name 'gms2.pl' -print -quit)
    if [[ -z "\${gms2_script}" ]]; then
        echo "gms2.pl was not found in the GeneMarkS-2 home" >&2
        exit 1
    fi
    gms2_script=\$(readlink -f "\${gms2_script}")
    gms2_dir=\$(dirname "\${gms2_script}")
    chmod u+x "\${gms2_script}" "\${gms2_dir}"/* 2>/dev/null || true

    install -m 600 "${genemark_key}" genemark_user/.gmhmmp2_key
    cp genemark_user/.gmhmmp2_key genemark_user/.gm_key

    GENEMARK_USER_DIR="\$PWD/genemark_user"
    export GMHMMP2_KEY="\${GENEMARK_USER_DIR}/.gmhmmp2_key"
    export PATH="\${gms2_dir}:\$PATH"
    export LC_ALL=C
    export LANG=C

    mag_path=\$(readlink -f "${mag}")
    output_root="\$PWD"

    echo "Running GeneMarkS-2 gene prediction"

    cd "\${GENEMARK_USER_DIR}"
    perl "\${gms2_script}" \
        --seq "\${mag_path}" \
        --genome-type bacteria \
        --threads ${task.cpus} \
        --format gff \
        --output "\${output_root}/${prefix}.raw.gff" \
        --fnn "\${output_root}/${prefix}.raw.fnn" \
        --faa "\${output_root}/${prefix}.raw.faa" \
        ${args} \
        > "\${output_root}/${prefix}.genemarks2.log" 2>&1
    cd "\${output_root}"

    awk -v MAG='${mag_id}' '
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
    ' "${prefix}.raw.faa" > "${prefix}.proteins.faa"

    awk -v MAG='${mag_id}' '
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
    ' "${prefix}.raw.fnn" > "${prefix}.cds.fna"

    awk -v MAG='${mag_id}' '
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
    ' "${prefix}.raw.gff" > "${prefix}.genes.gff"

    test -s "${prefix}.proteins.faa"
    test -s "${prefix}.cds.fna"
    test -s "${prefix}.genes.gff"
    """

    stub:
    def prefix = task.ext.prefix ?: (meta.mag_id ?: meta.id)
    def mag_id = meta.mag_id ?: meta.id

    """
    set -euo pipefail

    printf '>%s|stub_contig|1\nMKKIGYSAPRQTKEAIEKLA\n' \
        '${mag_id}' > "${prefix}.proteins.faa"
    printf '>%s|stub_contig|1\nATGAAAAAAATCGGTTATTCAGCTCCTCGTCAAACCAAAGAAGCTATTGAAAAACTGGCT\n' \
        '${mag_id}' > "${prefix}.cds.fna"
    printf '##gff-version 3\nstub_contig\tGeneMarkS-2\tCDS\t1\t63\t.\t+\t0\tID=%s|stub_contig|1;gene_id=%s|stub_contig|1;protein_id=%s|stub_contig|1;mag_id=%s\n' \
        '${mag_id}' '${mag_id}' '${mag_id}' '${mag_id}' > "${prefix}.genes.gff"

    printf 'GeneMarkS-2 stub gene prediction completed\n' > "${prefix}.genemarks2.log"
    """
}
