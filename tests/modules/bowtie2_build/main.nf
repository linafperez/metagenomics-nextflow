process BOWTIE2_BUILD_SYNTHETIC {
    tag "synthetic_host"
    label 'process_single'

    conda "${moduleDir}/../../../modules/core/bowtie2/environment.yml"
    container 'quay.io/biocontainers/bowtie2:2.5.4--he20e202_2'

    input:
    path reference

    output:
    path 'synthetic_host.*.bt2*', arity: '6', emit: index
    val 'synthetic_host', emit: prefix
    path 'synthetic_host.bowtie2-build.log', emit: log
    tuple val("${task.process}"), val('bowtie2'), val('2.5.4'), emit: versions

    script:
    """
    echo "Building the synthetic host Bowtie2 index"

    bowtie2-build \
        --threads ${task.cpus} \
        "${reference}" \
        synthetic_host \
        > synthetic_host.bowtie2-build.log 2>&1
    """

    stub:
    """
    for suffix in 1 2 3 4 rev.1 rev.2; do
        printf 'Bowtie2 2.5.4 synthetic index stub: %s\n' "\${suffix}" \
            > "synthetic_host.\${suffix}.bt2"
    done
    printf 'Bowtie2 synthetic index stub completed\n' \
        > synthetic_host.bowtie2-build.log
    """
}
