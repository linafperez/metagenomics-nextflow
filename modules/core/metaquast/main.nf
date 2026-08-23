process METAQUAST {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/quast:5.3.0--py313pl5321h5ca1c30_2'

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path('*.metaquast'), emit: results
    tuple val(meta), path('*.metaquast.report.tsv'), emit: report_tsv
    tuple val(meta), path('*.metaquast.report.html'), emit: report_html
    tuple val(meta), path('*.metaquast.log'), emit: log
    tuple val("${task.process}"), val('metaquast'), val('5.3.0'), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id

    """
    echo "Running MetaQUAST assembly evaluation"

    export MPLCONFIGDIR="\$PWD/.matplotlib"
    mkdir -p "\$MPLCONFIGDIR"

    metaquast.py \
        "${assembly}" \
        --output-dir "${prefix}.metaquast" \
        --threads ${task.cpus} \
        --max-ref-number 0 \
        ${args} \
        > "${prefix}.metaquast.log" 2>&1

    cp "${prefix}.metaquast/report.tsv" "${prefix}.metaquast.report.tsv"
    cp "${prefix}.metaquast/report.html" "${prefix}.metaquast.report.html"
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id

    """
    mkdir -p "${prefix}.metaquast"

    printf '%s\\t%s\\n' \
        'Assembly' '${prefix}.contigs.fa' \
        '# contigs (>= 0 bp)' '1' \
        '# contigs (>= 1000 bp)' '1' \
        'Total length (>= 0 bp)' '2432' \
        'Total length (>= 1000 bp)' '2432' \
        '# contigs' '1' \
        'Largest contig' '2432' \
        'Total length' '2432' \
        'GC (%)' '50.00' \
        'N50' '2432' \
        'N90' '2432' \
        'L50' '1' \
        'L90' '1' \
        > "${prefix}.metaquast/report.tsv"

    printf '%s\\n' \
        '<!doctype html>' \
        '<html lang="en"><head><meta charset="utf-8"><title>MetaQUAST stub report</title></head>' \
        '<body><h1>MetaQUAST stub report</h1><p>One synthetic contig, 2432 bp.</p></body></html>' \
        > "${prefix}.metaquast/report.html"

    cp "${prefix}.metaquast/report.tsv" "${prefix}.metaquast.report.tsv"
    cp "${prefix}.metaquast/report.html" "${prefix}.metaquast.report.html"
    printf 'MetaQUAST stub evaluation completed\n' > "${prefix}.metaquast.log"
    """
}
