process SRATOOLS_ACQUIRE {
    tag "${meta.id}:${run_accession}"
    label 'process_high'

    container "${params.sraContainer}"
    containerOptions "${params.sraContainerOptions ?: ''}"
    conda "${moduleDir}/environment.yml"
    scratch params.sraScratchDir ?: false

    input:
    tuple val(meta), val(run_accession)
    path pair_validator

    output:
    tuple val(meta), val(run_accession), path("${run_accession}_R*.fastq.gz", arity: 2), emit: reads
    tuple val(meta), path("${run_accession}.sra.log"), emit: log
    tuple val("${task.process}"), val('sra-tools'), val('3.4.1'), emit: versions

    script:
    """
    prefetch --version 2>&1 | grep -F '3.4.1' >/dev/null
    vdb-validate --version 2>&1 | grep -F '3.4.1' >/dev/null
    fasterq-dump --version 2>&1 | grep -F '3.4.1' >/dev/null
    pigz --version 2>&1 | grep -F '2.8' >/dev/null

    # All temporary directories belong exclusively to this task.
    mkdir -p sra_cache fasterq_tmp fastq task_home
    export HOME="\$PWD/task_home"
    export TMPDIR="\$PWD/fasterq_tmp"
    {
        prefetch '${run_accession}' --max-size '${params.sraMaxSize}' --output-directory sra_cache
        test -s 'sra_cache/${run_accession}/${run_accession}.sra'
        vdb-validate 'sra_cache/${run_accession}'
        fasterq-dump 'sra_cache/${run_accession}' --split-files --skip-technical \\
            --threads ${task.cpus} --temp fasterq_tmp --outdir fastq

        test -s 'fastq/${run_accession}_1.fastq'
        test -s 'fastq/${run_accession}_2.fastq'
        if [[ -e 'fastq/${run_accession}.fastq' ]]; then
            echo 'Discarding singleton reads; paired-end analysis retains mates only.' >&2
            rm -- 'fastq/${run_accession}.fastq'
        fi
        awk -v mate2='fastq/${run_accession}_2.fastq' -f '${pair_validator}' 'fastq/${run_accession}_1.fastq'
        pigz -f -p ${task.cpus} 'fastq/${run_accession}_1.fastq'
        pigz -f -p ${task.cpus} 'fastq/${run_accession}_2.fastq'
        pigz -t 'fastq/${run_accession}_1.fastq.gz' 'fastq/${run_accession}_2.fastq.gz'
        mv -- 'fastq/${run_accession}_1.fastq.gz' '${run_accession}_R1.fastq.gz'
        mv -- 'fastq/${run_accession}_2.fastq.gz' '${run_accession}_R2.fastq.gz'
    } > '${run_accession}.sra.log' 2>&1

    # Clean task-local temporaries after validation; keep cached outputs in work/.
    rm -rf -- sra_cache fasterq_tmp fastq task_home
    """
}
