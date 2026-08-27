process SEMIBIN2 {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/semibin:1.5.0--pyhdfd78af_1'

    input:
    tuple val(meta), path(contigs), path(bams)

    output:
    tuple val(meta), path('*.bins/*.{fa,fna,fasta}', arity: '1..*'), emit: bins
    tuple val(meta), path('*.contigs2bin.tsv'), emit: contigs2bin
    tuple val(meta), path('*.semibin2'), emit: native_outputs
    tuple val(meta), path('*.semibin2.log'), emit: log
    tuple val("${task.process}"), val('semibin'), val('1.5.0'), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: meta.id
    def bam_files = bams instanceof List ? bams : [bams]
    if (!bam_files) {
        error 'SEMIBIN2 requires at least one reference-sorted BAM file'
    }
    def bam_args = bam_files.collect { bam -> "\"${bam}\"" }.join(' ')
    def standardize_bins = '''
    shopt -s nullglob
    native_bins=("$NATIVE_BINS"/*.fa "$NATIVE_BINS"/*.fna "$NATIVE_BINS"/*.fasta)
    if [ "${#native_bins[@]}" -eq 0 ]; then
        echo "SemiBin2 did not produce any FASTA bins" >&2
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

    mkdir -p "${prefix}.semibin2.bins"
    SemiBin2 single_easy_bin \
        -i "${contigs}" \
        -b ${bam_args} \
        -t ${task.cpus} \
        --engine cpu \
        --output-compression none \
        -o "${prefix}.semibin2" \
        ${args} \
        2> >(tee "${prefix}.semibin2.log" >&2)

    NATIVE_BINS="${prefix}.semibin2/output_bins"
    BINS_DIR="${prefix}.semibin2.bins"
    MAP_FILE="${prefix}.semibin2.contigs2bin.tsv"
    ${standardize_bins}

    """

}
