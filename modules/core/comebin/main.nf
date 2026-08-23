process COMEBIN {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/comebin:1.0.4--hdfd78af_1'

    input:
    tuple val(meta), path(contigs), path(bams)

    output:
    tuple val(meta), path('*.bins/*.{fa,fna,fasta}', arity: '1..*'), emit: bins
    tuple val(meta), path('*.contigs2bin.tsv'), emit: contigs2bin
    tuple val(meta), path('*.comebin'), emit: native_outputs
    tuple val(meta), path('*.comebin.log'), emit: log
    tuple val("${task.process}"), val('comebin'), val('1.0.4'), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: meta.id
    def bam_files = bams instanceof List ? bams : [bams]
    if (!bam_files) {
        error 'COMEBIN requires at least one reference-sorted BAM file'
    }
    def bam_links = bam_files.collect { bam -> "ln -s \"\$(readlink -f '${bam}')\" \"bam/${bam.name}\"" }.join('\n')
    def standardize_bins = '''
    shopt -s nullglob
    native_bins=("$NATIVE_BINS"/*.fa "$NATIVE_BINS"/*.fna "$NATIVE_BINS"/*.fasta)
    if [ "${#native_bins[@]}" -eq 0 ]; then
        echo "COMEBin did not produce any FASTA bins" >&2
        exit 1
    fi
    cp -- "${native_bins[@]}" "$BINS_DIR/"

    : > "$MAP_FILE"
    bin_files=("$BINS_DIR"/*.fa "$BINS_DIR"/*.fna "$BINS_DIR"/*.fasta)
    for bin_file in "${bin_files[@]}"; do
        bin_name=$(basename "$bin_file")
        bin_name=${bin_name%.*}
        awk -v bin="$bin_name" 'BEGIN { OFS="\t" } /^>/ { sub(/^>/, ""); split($0, fields, /[[:space:]]+/); print fields[1], bin }' "$bin_file" >> "$MAP_FILE"
    done
    test -s "$MAP_FILE"
    '''.stripIndent()

    """
    set -euo pipefail

    mkdir -p bam "${prefix}.comebin.bins"
    ${bam_links}

    CUDA_VISIBLE_DEVICES='' run_comebin.sh \
        -a "${contigs}" \
        -o "${prefix}.comebin" \
        -p bam \
        -t ${task.cpus} \
        ${args} \
        2> >(tee "${prefix}.comebin.log" >&2)

    NATIVE_BINS="${prefix}.comebin/comebin_res/comebin_res_bins"
    BINS_DIR="${prefix}.comebin.bins"
    MAP_FILE="${prefix}.comebin.contigs2bin.tsv"
    ${standardize_bins}

    """

    stub:
    def prefix = task.ext.prefix ?: meta.id

    """
    set -euo pipefail

    mkdir -p "${prefix}.comebin.bins" "${prefix}.comebin/comebin_res/comebin_res_bins"
    awk 'BEGIN { found=0 } /^>/ { if (found) exit; found=1 } found { print }' "${contigs}" > "${prefix}.comebin.bins/${prefix}.comebin.1.fa"
    cp "${prefix}.comebin.bins/${prefix}.comebin.1.fa" "${prefix}.comebin/comebin_res/comebin_res_bins/"
    awk -v bin="${prefix}.comebin.1" 'BEGIN { OFS="\t" } /^>/ { sub(/^>/, ""); split(\$0, fields, /[[:space:]]+/); print fields[1], bin; exit }' "${contigs}" > "${prefix}.comebin.contigs2bin.tsv"
    cp "${prefix}.comebin.contigs2bin.tsv" "${prefix}.comebin/comebin_res/comebin_res.tsv"
    printf 'COMEBin stub completed\n' > "${prefix}.comebin.log"

    """
}
