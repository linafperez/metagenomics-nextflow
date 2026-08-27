process METABAT2 {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/metabat2:2.18--h38e344b_2'

    input:
    tuple val(meta), path(contigs), path(metabat_depth)

    output:
    tuple val(meta), path('*.bins/*.{fa,fna,fasta}', arity: '1..*'), emit: bins
    tuple val(meta), path('*.contigs2bin.tsv'), emit: contigs2bin
    tuple val(meta), path('*.metabat2.log'), emit: log
    tuple val("${task.process}"), val('metabat2'), val('2.18'), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id
    def make_map = '''
    shopt -s nullglob
    bin_files=("$BINS_DIR"/*.fa "$BINS_DIR"/*.fna "$BINS_DIR"/*.fasta)
    if [ "${#bin_files[@]}" -eq 0 ]; then
        echo "MetaBAT2 did not produce any FASTA bins" >&2
        exit 1
    fi

    : > "$MAP_FILE"
    for bin_file in "${bin_files[@]}"; do
        bin_name=$(basename "$bin_file")
        bin_name=${bin_name%.*}
        awk -v bin="$bin_name" 'BEGIN { OFS="\t" } /^>/ { sub(/^>/, ""); split($0, fields, /[[:space:]]+/); print fields[1], bin }' "$bin_file" >> "$MAP_FILE"
    done
    test -s "$MAP_FILE"
    '''.stripIndent()

    """
    set -euo pipefail

    mkdir -p "${prefix}.metabat2.bins"
    metabat2 \
        -i "${contigs}" \
        -a "${metabat_depth}" \
        -o "${prefix}.metabat2.bins/${prefix}.metabat2" \
        -t ${task.cpus} \
        ${args} \
        2> >(tee "${prefix}.metabat2.log" >&2)

    BINS_DIR="${prefix}.metabat2.bins"
    MAP_FILE="${prefix}.metabat2.contigs2bin.tsv"
    ${make_map}

    """

}
