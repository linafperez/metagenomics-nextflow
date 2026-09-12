process SRATOOLS_ACQUIRE {
    tag "${meta.id}"
    label 'process_high'

    container "${params.sraContainer}"
    containerOptions "${params.sraContainerOptions ?: ''}"
    conda "${moduleDir}/environment.yml"

    input:
    val meta
    path run_manifest
    path acquisition_helper
    val scratch_root
    val cache_root
    val temporary_root
    val maximum_size

    output:
    tuple val(meta), path("${meta.id}_R*.fastq.gz", arity: 2), emit: reads
    path 'sra_acquisition_versions.tsv', emit: versions

    script:
    """
    prefetch --version 2>&1 | grep -F '3.4.1' >/dev/null
    vdb-validate --version 2>&1 | grep -F '3.4.1' >/dev/null
    fasterq-dump --version 2>&1 | grep -F '3.4.1' >/dev/null
    pigz --version 2>&1 | grep -F '2.8' >/dev/null

    python3 "${acquisition_helper}" \
        --manifest "${run_manifest}" \
        --sample-id "${meta.id}" \
        --output-dir . \
        --scratch-dir "${scratch_root}" \
        --prefetch-dir "${cache_root}" \
        --temp-dir "${temporary_root}" \
        --threads ${task.cpus} \
        --max-size "${maximum_size}" \
        --force

    printf 'SRATOOLS_ACQUIRE\tsra-tools\t3.4.1\nSRATOOLS_ACQUIRE\tpigz\t2.8\nSRATOOLS_ACQUIRE\tpython\t3.12.11\n' \
        > sra_acquisition_versions.tsv
    """
}
