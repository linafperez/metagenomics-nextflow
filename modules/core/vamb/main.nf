process VAMB {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/vamb:5.0.4--pyhdfd78af_0'

    input:
    tuple val(meta), path(contigs), path(vamb_abundance)

    output:
    tuple val(meta), path('*.bins/*.{fa,fna,fasta}', arity: '1..*'), emit: bins
    tuple val(meta), path('*.contigs2bin.tsv'), emit: contigs2bin
    tuple val(meta), path('*.vamb'), emit: native_outputs
    tuple val(meta), path('*.vamb.log'), emit: log
    tuple val("${task.process}"), val('vamb'), val('5.0.4'), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id
    def standardize_bins = '''
    shopt -s nullglob
    native_bins=("$NATIVE_BINS"/*.fa "$NATIVE_BINS"/*.fna "$NATIVE_BINS"/*.fasta)
    if [ "${#native_bins[@]}" -eq 0 ]; then
        echo "Vamb did not materialize any FASTA bins" >&2
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

    mkdir -p "${prefix}.vamb.bins"
    vamb bin default \
        --fasta "${contigs}" \
        --abundance_tsv "${vamb_abundance}" \
        --outdir "${prefix}.vamb" \
        -m 2000 \
        --minfasta 1 \
        -e 100 \
        -q 25 75 \
        -p ${task.cpus} \
        ${args} \
        2> >(tee "${prefix}.vamb.log" >&2)

    NATIVE_BINS="${prefix}.vamb/bins"
    BINS_DIR="${prefix}.vamb.bins"
    MAP_FILE="${prefix}.vamb.contigs2bin.tsv"
    ${standardize_bins}
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id

    """
    set -euo pipefail

    mkdir -p "${prefix}.vamb.bins" "${prefix}.vamb/bins"
    awk 'BEGIN { found=0 } /^>/ { if (found) exit; found=1 } found { print }' "${contigs}" > "${prefix}.vamb.bins/${prefix}.vamb.1.fna"
    cp "${prefix}.vamb.bins/${prefix}.vamb.1.fna" "${prefix}.vamb/bins/"
    contig_id=\$(awk '/^>/{sub(/^>/, ""); split(\$0, fields, /[[:space:]]+/); print fields[1]; exit}' "${contigs}")
    contig_id=\${contig_id:-stub_contig}
    printf '%s\t%s\n' "\${contig_id}" "${prefix}.vamb.1" > "${prefix}.vamb.contigs2bin.tsv"
    printf 'clustername\tcontigname\n%s\t%s\n' "${prefix}.vamb.1" "\${contig_id}" > "${prefix}.vamb/vae_clusters_unsplit.tsv"
    printf 'Name\tRadius\tPeak/valley ratio\tKind\tBp\tNcontigs\tMedoid\n%s\t0\t0\tNormal\t1000\t1\t%s\n' "${prefix}.vamb.1" "\${contig_id}" > "${prefix}.vamb/vae_clusters_metadata.tsv"
    printf 'Vamb stub completed\n' > "${prefix}.vamb/log.txt"
    cp "${prefix}.vamb/log.txt" "${prefix}.vamb.log"
    """
}
