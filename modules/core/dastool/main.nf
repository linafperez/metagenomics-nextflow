process DASTOOL {
    tag "${meta.id}"
    label 'process_high_memory'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/das_tool:1.1.7--r44hdfd78af_1'

    input:
    tuple val(meta), path(contigs), path(comebin_map), path(metabat2_map), path(semibin2_map), path(vamb_map)

    output:
    tuple val(meta), path('*_DASTool_bins/*.{fa,fna,fasta}', arity: '1..*'), emit: bins
    tuple val(meta), path('*_DASTool_summary.tsv'), emit: summary
    tuple val(meta), path('*_DASTool_contig*bin.tsv'), emit: contigs2bin
    tuple val(meta), path('*_allBins.eval'), emit: evaluations
    tuple val(meta), path('*_DASTool.log'), emit: log
    tuple val(meta), path('*_proteins.faa'), optional: true, emit: proteins
    tuple val("${task.process}"), val('das_tool'), val('1.1.7'), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id
    def keep_native_outputs = params.save_intermediates.toString().toBoolean()

    """
    set -euo pipefail

    for bin_map in "${comebin_map}" "${metabat2_map}" "${semibin2_map}" "${vamb_map}"; do
        if [ ! -s "\${bin_map}" ]; then
            echo "DAS Tool input map is empty: \${bin_map}" >&2
            exit 1
        fi
        awk 'NF < 2 { exit 1 }' "\${bin_map}" || {
            echo "DAS Tool input map is not a two-column table: \${bin_map}" >&2
            exit 1
        }
    done

    DAS_Tool \
        --bins "${comebin_map},${metabat2_map},${semibin2_map},${vamb_map}" \
        --labels COMEBin,MetaBAT2,SemiBin2,Vamb \
        --contigs "${contigs}" \
        --outputbasename "${prefix}" \
        --threads ${task.cpus} \
        --search_engine diamond \
        --write_bins \
        --write_bin_evals \
        ${args}

    test -n "\$(find "${prefix}_DASTool_bins" -type f \( -name '*.fa' -o -name '*.fna' -o -name '*.fasta' \) -size +0c -print -quit)"
    test -s "${prefix}_DASTool_summary.tsv"
    test -s "${prefix}_DASTool_contigs2bin.tsv"
    test -s "${prefix}_allBins.eval"
    test -s "${prefix}_DASTool.log"

    rm -f -- \
        "${prefix}.seqlength" \
        "${prefix}_proteins.faa.bacteria.scg" \
        "${prefix}_proteins.faa.archaea.scg" \
        "${prefix}_proteins.faa.findSCG.b6" \
        "${prefix}_proteins.faa.all.b6" \
        "${prefix}_proteins.faa.scg.candidates.faa" \
        all_prot.dmnd
    ${keep_native_outputs ? '' : "rm -f -- '${prefix}_proteins.faa'"}
    """

}
