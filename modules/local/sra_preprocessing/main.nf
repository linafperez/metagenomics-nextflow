process SRA_ACQUIRE {
    tag "${sample_id}"
    label 'process_high'

    container "${params.sraContainer}"
    containerOptions "${params.sraContainerOptions ?: ''}"
    conda "${moduleDir}/environment.yml"

    input:
    val sample_id
    path run_manifest
    path acquisition_helper
    val scratch_root
    val cache_root
    val temporary_root
    val maximum_size

    output:
    tuple val(sample_id), path("${sample_id}_R*.fastq.gz", arity: 2), emit: reads
    path 'sra_acquisition_versions.tsv', emit: versions

    script:
    """
    prefetch --version 2>&1 | grep -F '3.4.1' >/dev/null
    fasterq-dump --version 2>&1 | grep -F '3.4.1' >/dev/null
    pigz --version 2>&1 | grep -F '2.8' >/dev/null

    python3 "${acquisition_helper}" \
        --manifest "${run_manifest}" \
        --sample-id "${sample_id}" \
        --output-dir . \
        --scratch-dir "${scratch_root}" \
        --prefetch-dir "${cache_root}" \
        --temp-dir "${temporary_root}" \
        --threads ${task.cpus} \
        --max-size "${maximum_size}" \
        --force

    printf 'SRA_ACQUIRE\tsra-tools\t3.4.1\nSRA_ACQUIRE\tpigz\t2.8\nSRA_ACQUIRE\tpython\t3.12.11\n' \
        > sra_acquisition_versions.tsv
    """
}

process PERSIST_SRA_CHECKPOINT {
    tag "${meta.id}"
    label 'process_medium'
    // The external checkpoint root is mutable state outside the task hash.
    // Always replay this durable commit instead of trusting a cached side effect.
    cache false

    container 'python:3.12.11-slim-bookworm'
    containerOptions "${params.sraContainerOptions ?: ''}"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(reads, arity: 2)
    path run_manifest
    path reports
    val version_records
    path checkpoint_helper
    val checkpoint_root

    output:
    path 'checkpoint_record.json', emit: record
    path 'preprocessing_versions.tsv', emit: versions

    script:
    def version_text = version_records
        .collect { record -> record.collect { value -> value.toString() }.join('\t') }
        .sort()
        .join('\n')
    def report_arguments = reports
        .collect { report -> "--report '${report}'" }
        .join(' ')

    """
    cat > preprocessing_versions.tsv <<'NF_SRA_VERSIONS'
${version_text}
NF_SRA_VERSIONS

    python3 "${checkpoint_helper}" persist \
        --run-manifest "${run_manifest}" \
        --sample-id "${meta.id}" \
        --read-1 "${reads[0]}" \
        --read-2 "${reads[1]}" \
        --checkpoint-dir "${checkpoint_root}" \
        ${report_arguments} \
        --report preprocessing_versions.tsv \
        --output-record checkpoint_record.json
    """
}

process CHECK_SRA_CHECKPOINTS {
    tag 'frozen SRA cohort'
    label 'process_single'
    // Always observe current external checkpoint contents on staged resumes.
    cache false

    container 'python:3.12.11-slim-bookworm'
    containerOptions "${params.sraContainerOptions ?: ''}"
    conda "${moduleDir}/environment.yml"

    input:
    path run_manifest
    path checkpoint_helper
    val checkpoint_root
    val require_complete

    output:
    path 'sra_checkpoint_manifest.tsv', emit: manifest
    path 'sra_pending_samples.tsv', emit: pending
    path 'sra_checkpoint_status.json', emit: status

    script:
    def complete_argument = require_complete ? '--require-complete' : ''
    """
    python3 "${checkpoint_helper}" reconcile \
        --run-manifest "${run_manifest}" \
        --checkpoint-dir "${checkpoint_root}" \
        --output-manifest sra_checkpoint_manifest.tsv \
        --pending-output sra_pending_samples.tsv \
        --status-output sra_checkpoint_status.json \
        ${complete_argument}
    """
}

process FINALIZE_SRA_GLOBAL_RUN {
    tag "${project_accession}"
    label 'process_single'

    container 'python:3.12.11-slim-bookworm'
    conda "${moduleDir}/environment.yml"

    input:
    val project_accession
    path checkpoint_manifest
    path multiqc_report
    path software_versions
    path mag_abundance
    val durable_checkpoint_manifest_path
    val durable_multiqc_path
    val durable_versions_path
    val durable_abundance_path

    output:
    path 'sra_global_success.json', emit: marker

    script:
    """
    test -s "${multiqc_report}"
    test -s "${software_versions}"
    test -s "${mag_abundance}"
    python3 - <<'PY'
import json
import hashlib
import os
from datetime import datetime, timezone

def describe(source, durable_path):
    digest = hashlib.sha256()
    with open(source, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return {
        'path': durable_path,
        'bytes': os.path.getsize(source),
        'sha256': digest.hexdigest(),
    }

with open('sra_global_success.json', 'w', encoding='utf-8') as handle:
    json.dump({
        'schema_version': 1,
        'project_accession': '${project_accession}',
        'status': 'complete',
        'completed_at_utc': datetime.now(timezone.utc).isoformat(timespec='seconds').replace('+00:00', 'Z'),
        'checkpoint_manifest': describe('${checkpoint_manifest}', '${durable_checkpoint_manifest_path}'),
        'outputs': {
            'multiqc_report': describe('${multiqc_report}', '${durable_multiqc_path}'),
            'software_versions': describe('${software_versions}', '${durable_versions_path}'),
            'mag_abundance': describe('${mag_abundance}', '${durable_abundance_path}'),
        },
    }, handle, indent=2, sort_keys=True)
    handle.write('\n')
PY
    """
}
