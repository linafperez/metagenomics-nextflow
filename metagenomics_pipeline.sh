#!/usr/bin/env bash

set -euo pipefail

readonly MINIMUM_NEXTFLOW_VERSION="26.04.6"
readonly MINIMUM_PYTHON_VERSION="3.10.0"
readonly PIPELINE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
readonly DATABASE_PREPARER="${PIPELINE_ROOT}/bin/prepare_databases.sh"
readonly CHECKPOINT_MANAGER="${PIPELINE_ROOT}/bin/manage_sra_checkpoints.py"
readonly SRA_RESOLVER="${PIPELINE_ROOT}/bin/resolve_sra_project.py"
readonly STORAGE_MONITOR="${PIPELINE_ROOT}/bin/monitor_storage.py"
readonly RESOURCE_SUMMARIZER="${PIPELINE_ROOT}/bin/summarize_resources.py"
readonly SLURM_ACCOUNTING_COLLECTOR="${PIPELINE_ROOT}/bin/collect_slurm_accounting.py"

environment=""
runtime=""
mode=""
database_config=""
db_root=""
dry_run=false
resume=false
storage_constrained=false
enable_gpu=false
input_path=""
sra_project=""
sra_checkpoint_dir=""
sra_scratch_root=""
sra_cache_dir=""
sra_temp_dir=""
sra_email=""
sra_platforms="ILLUMINA,BGISEQ"
sra_max_size="u"
keep_sra_checkpoints=false
outdir="results"
work_dir=""
resource_database_root=""
resource_sample_interval="60"
gpu_accelerators="1"
gpu_telemetry_interval="10"
slurm_gpu_gres=""
monitor_pid=""
monitor_stop_file=""
monitor_sample_request_file=""
monitor_stage_file=""
telemetry_started=false
project_status="not_started"
project_exit_code=0
project_manifest=""
trace_registry=""
project_status_file=""
resource_root=""
run_token=""
resume_registry=""
publish_dir_mode_override=""
run_lock_token=""
declare -a forwarded_args=()
declare -a profiles=()
declare -a run_lock_dirs=()

usage() {
    cat <<'USAGE'
Usage:
  ./metagenomics_pipeline.sh --<environment> --<runtime> --run --input FILE [options]
  ./metagenomics_pipeline.sh --<environment> --<runtime> --run \
      --sra-project PRJNA... --sra-checkpoint-dir PATH --sra-scratch-dir PATH [options]
  ./metagenomics_pipeline.sh --<environment> --prepare-databases --db-root PATH [options]

Environment (select one):
  --local                  Use the local executor.
  --hpc                    Use the SLURM executor.

Software runtime (select one for pipeline execution):
  --docker | --conda | --apptainer | --singularity

Mode (select one):
  --run                    Run a production workflow.
  --prepare-databases      Prepare shared databases separately.

Input (select exactly one for --run):
  --input FILE             Existing paired-FASTQ samplesheet.
  --sra-project ACCESSION  Frozen BioProject mode (PRJNA/PRJEB/PRJDB).

SRA lifecycle options:
  --sra-checkpoint-dir PATH  Required durable non-host-read checkpoint root.
  --sra-scratch-dir PATH     Required disposable SRA/Nextflow scratch root.
  --sra-cache-dir PATH       Prefetch cache root (default: SCRATCH/sra-cache).
  --sra-temp-dir PATH        fasterq temporary root (default: SCRATCH/fasterq-temp).
  --sra-email EMAIL          Optional NCBI E-utilities contact email.
  --sra-platforms LIST       Allowed platforms (default: ILLUMINA,BGISEQ).
  --sra-max-size VALUE       prefetch maximum per run (default: u).
  --keep-sra-checkpoints     Retain non-host FASTQ checkpoints after global success.

Execution and accounting:
  --database-config FILE
  --work-dir PATH             Work root for local FASTQ mode.
  --outdir PATH               Results root (default: results).
  --storage-constrained       Serialize local-FASTQ tasks; SRA mode always does this.
  --enable-gpu                Enable verified COMEBin/SemiBin2/Vamb GPU paths.
  --gpu-accelerators N        GPUs requested by each enabled process (default: 1).
  --gpu-telemetry-interval N   GPU sampling interval in seconds (default: 10).
  --slurm-gpu-gres VALUE      Site-specific SLURM GRES value, e.g. gpu:a100:1.
  --resource-sample-interval SECONDS
  --resource-database-root PATH
  --resume
  --dry-run
  -h, --help
  --version

Unrecognized arguments are forwarded unchanged to every applicable Nextflow
stage. BioProject mode never downloads the entire project at once: the launcher
waits for a durable checkpoint before starting the next biological sample.
USAGE
}

die() {
    printf 'Error: %s\n' "$*" >&2
    exit 2
}

set_selection() {
    local selection_name="$1" selection_value="$2" current_value="$3"
    if [[ -n "${current_value}" && "${current_value}" != "${selection_value}" ]]; then
        die "select exactly one ${selection_name}"
    fi
}

require_value() {
    local option="$1" value="${2:-}"
    [[ -n "${value}" && "${value}" != --* ]] || die "${option} requires a value"
}

print_command() {
    local argument
    printf 'Command:'
    for argument in "$@"; do
        printf ' %q' "${argument}"
    done
    printf '\n'
}

absolute_path() {
    local value="$1"
    if command -v realpath >/dev/null 2>&1; then
        realpath -m -- "${value}"
    elif [[ "${value}" == /* ]]; then
        printf '%s\n' "${value}"
    else
        printf '%s/%s\n' "${PWD}" "${value}"
    fi
}

path_is_within() {
    local candidate="$1" parent="$2"
    [[ "${candidate}" == "${parent}" || "${candidate}" == "${parent}/"* ]]
}

require_container_safe_path() {
    local option="$1" path="$2"
    [[ ! "${path}" =~ [[:space:]] && "${path}" != *","* && "${path}" != *":"* ]] \
        || die "${option} cannot contain whitespace, comma, or colon in container-backed SRA mode"
}

require_interpolation_safe_path() {
    local option="$1" path="$2"
    [[ "${path}" != *"'"* && "${path}" != *'"'* && "${path}" != *'`'* \
        && "${path}" != *'$'* && "${path}" != *'\'* \
        && "${path}" != *$'\n'* && "${path}" != *$'\r'* ]] \
        || die "${option} contains characters unsafe for staged script interpolation"
}

release_run_locks() {
    local index lock_dir claim_file metadata_file metadata_temporary claim_value
    for ((index = ${#run_lock_dirs[@]} - 1; index >= 0; index--)); do
        lock_dir="${run_lock_dirs[index]}"
        claim_file="${lock_dir}/.owner-${run_lock_token}"
        metadata_file="${lock_dir}/owner.tsv"
        metadata_temporary="${lock_dir}/.owner.tsv.${run_lock_token}.tmp"

        if [[ -L "${lock_dir}" || ! -d "${lock_dir}" ]]; then
            printf 'Warning: run lock path changed before release; it was left untouched: %s\n' \
                "${lock_dir}" >&2
            continue
        fi
        claim_value=""
        if [[ ! -f "${claim_file}" || -L "${claim_file}" ]] \
            || ! IFS= read -r claim_value < "${claim_file}" \
            || [[ "${claim_value}" != "${run_lock_token}" ]]; then
            printf 'Warning: run lock ownership marker is missing or changed; lock was left untouched: %s\n' \
                "${lock_dir}" >&2
            continue
        fi
        if ! rm -f -- "${metadata_temporary}" "${metadata_file}" "${claim_file}"; then
            printf 'Warning: run lock metadata could not be removed; lock was left for inspection: %s\n' \
                "${lock_dir}" >&2
            continue
        fi
        if ! rmdir -- "${lock_dir}" 2>/dev/null; then
            printf 'Warning: run lock contains unexpected content and was left for inspection: %s\n' \
                "${lock_dir}" >&2
        fi
    done
    run_lock_dirs=()
}

acquire_run_lock() {
    local lock_dir="$1" scope="$2" checkpoint_root="${3:-}"
    local claim_file metadata_file metadata_temporary lock_hostname

    if ! mkdir -- "${lock_dir}" 2>/dev/null; then
        if [[ -e "${lock_dir}" || -L "${lock_dir}" ]]; then
            die "run lock already exists at ${lock_dir}; another launcher may be active. This is fail-closed: inspect owner.tsv and remove the lock only after confirming that no run is active"
        fi
        die "could not create run lock at ${lock_dir}"
    fi
    claim_file="${lock_dir}/.owner-${run_lock_token}"
    metadata_file="${lock_dir}/owner.tsv"
    metadata_temporary="${lock_dir}/.owner.tsv.${run_lock_token}.tmp"
    if ! printf '%s\n' "${run_lock_token}" > "${claim_file}"; then
        rm -f -- "${claim_file}" 2>/dev/null || true
        rmdir -- "${lock_dir}" 2>/dev/null || true
        die "could not establish ownership of run lock ${lock_dir}"
    fi
    run_lock_dirs+=("${lock_dir}")
    lock_hostname="$(hostname 2>/dev/null || uname -n 2>/dev/null || printf 'unknown')"
    {
        printf 'field\tvalue\n'
        printf 'lock_token\t%s\n' "${run_lock_token}"
        printf 'scope\t%s\n' "${scope}"
        printf 'pid\t%s\n' "$$"
        printf 'hostname\t%q\n' "${lock_hostname}"
        printf 'started_at_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
        printf 'results_root\t%q\n' "${outdir}"
        printf 'checkpoint_root\t%q\n' "${checkpoint_root}"
    } > "${metadata_temporary}" || die "could not write run lock metadata at ${lock_dir}"
    mv -- "${metadata_temporary}" "${metadata_file}" \
        || die "could not publish run lock metadata at ${lock_dir}"
}

initialize_project_telemetry() {
    mkdir -p -- "${resource_root}/resources" "${resource_root}/invocations"
    if [[ ! -s "${trace_registry}" ]]; then
        printf 'invocation_id\ttrace_path\tsession_id\tstage\tstarted_at\tfinished_at\ttrace_raw\tlaunch_dir\tstatus\texit_code\n' > "${trace_registry}"
    fi
    if [[ ! -s "${resume_registry}" ]]; then
        printf 'resume_key\tsession_id\trecorded_at_utc\n' > "${resume_registry}"
    fi
    project_status="running"
    write_project_status "${project_status}" 0
    trap 'finish_telemetry $?' EXIT
    trap 'exit 130' INT TERM
}

version_at_least() {
    local actual="$1" required="$2" actual_part required_part index
    local -a actual_parts required_parts
    IFS=. read -r -a actual_parts <<< "${actual}"
    IFS=. read -r -a required_parts <<< "${required}"
    for index in 0 1 2; do
        actual_part="${actual_parts[index]:-0}"
        required_part="${required_parts[index]:-0}"
        actual_part="${actual_part%%[^0-9]*}"
        required_part="${required_part%%[^0-9]*}"
        actual_part="${actual_part:-0}"
        required_part="${required_part:-0}"
        ((10#${actual_part} > 10#${required_part})) && return 0
        ((10#${actual_part} < 10#${required_part})) && return 1
    done
    return 0
}

check_nextflow() {
    local detected_version
    command -v nextflow >/dev/null 2>&1 || die "Nextflow is not available on PATH"
    detected_version="$({ nextflow -version 2>&1 || true; } \
        | sed -nE 's/.*version[[:space:]]+([0-9]+\.[0-9]+\.[0-9]+).*/\1/p' \
        | head -n 1)"
    [[ -n "${detected_version}" ]] || die "could not determine the Nextflow version"
    version_at_least "${detected_version}" "${MINIMUM_NEXTFLOW_VERSION}" \
        || die "Nextflow ${MINIMUM_NEXTFLOW_VERSION} or newer is required; found ${detected_version}"
}

check_python() {
    local detected_version
    command -v python3 >/dev/null 2>&1 \
        || die "Python ${MINIMUM_PYTHON_VERSION} or newer is required for launcher helpers"
    detected_version="$(python3 -c 'import sys; print(".".join(map(str, sys.version_info[:3])))')" \
        || die "could not determine the Python version"
    version_at_least "${detected_version}" "${MINIMUM_PYTHON_VERSION}" \
        || die "Python ${MINIMUM_PYTHON_VERSION} or newer is required; found ${detected_version}"
}

check_runtime() {
    case "${runtime}" in
        docker)
            command -v docker >/dev/null 2>&1 || die "Docker is not available on PATH"
            docker info >/dev/null 2>&1 || die "the Docker daemon is not available"
            ;;
        conda)
            command -v mamba >/dev/null 2>&1 || die "Mamba is required by the Conda profile"
            ;;
        apptainer)
            command -v apptainer >/dev/null 2>&1 || die "Apptainer is not available on PATH"
            ;;
        singularity)
            command -v singularity >/dev/null 2>&1 || die "Singularity is not available on PATH"
            ;;
        *) die "unsupported software runtime '${runtime}'" ;;
    esac
}

write_project_status() {
    local status="$1" exit_code="$2" temporary="${project_status_file}.tmp.$$"
    [[ -n "${project_status_file}" ]] || return 0
    printf '{\n  "status": "%s",\n  "exit_code": %s,\n  "updated_at_utc": "%s"\n}\n' \
        "${status}" "${exit_code}" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "${temporary}"
    mv -f -- "${temporary}" "${project_status_file}"
}

ensure_trace_file() {
    local trace_file="$1"
    if [[ ! -s "${trace_file}" ]]; then
        printf 'task_id\thash\tnative_id\tprocess\tstatus\texit\tworkdir\n' > "${trace_file}"
    fi
}

extract_session_id() {
    local telemetry_dir="$1"
    [[ -s "${telemetry_dir}/nextflow.log" ]] || return 0
    sed -nE 's/.*Session UUID:[[:space:]]*([[:alnum:]-]+).*/\1/p' \
        "${telemetry_dir}/nextflow.log" | tail -n 1
}

latest_resume_session() {
    local resume_key="$1" session_id=""
    [[ -r "${resume_registry}" ]] || return 0
    session_id="$(awk -F '\t' -v key="${resume_key}" \
        'NR > 1 && $1 == key { session = $2 } END { print session }' \
        "${resume_registry}")"
    if [[ "${session_id}" =~ ^[0-9A-Fa-f]{8}-[0-9A-Fa-f]{4}-[0-9A-Fa-f]{4}-[0-9A-Fa-f]{4}-[0-9A-Fa-f]{12}$ ]]; then
        printf '%s\n' "${session_id}"
    elif [[ -n "${session_id}" ]]; then
        printf 'Warning: ignored invalid stored Nextflow session UUID for %s.\n' "${resume_key}" >&2
    fi
}

record_resume_session() {
    local resume_key="$1" telemetry_dir="$2" session_id
    session_id="$(extract_session_id "${telemetry_dir}")"
    if [[ "${session_id}" =~ ^[0-9A-Fa-f]{8}-[0-9A-Fa-f]{4}-[0-9A-Fa-f]{4}-[0-9A-Fa-f]{4}-[0-9A-Fa-f]{12}$ ]]; then
        printf '%s\t%s\t%s\n' "${resume_key}" "${session_id}" \
            "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "${resume_registry}"
    else
        printf 'Warning: no valid Nextflow session UUID was recorded for %s.\n' "${resume_key}" >&2
    fi
}

append_registry() {
    local invocation_id="$1" stage="$2" telemetry_dir="$3" started="$4" finished="$5" status="$6" exit_code="$7"
    local trace_file="${telemetry_dir}/execution_trace.tsv" session_id=""
    ensure_trace_file "${trace_file}"
    session_id="$(extract_session_id "${telemetry_dir}")"
    printf '%s\t%s\t%s\t%s\t%s\t%s\ttrue\t%s\t%s\t%s\n' \
        "${invocation_id}" "${trace_file}" "${session_id}" "${stage}" \
        "${started}" "${finished}" "${PIPELINE_ROOT}" "${status}" "${exit_code}" \
        >> "${trace_registry}"
}

run_nextflow() {
    local invocation_id="$1" stage="$2" invocation_work="$3"
    shift 3
    local telemetry_dir="${resource_root}/invocations/${run_token}_${invocation_id}"
    local started finished status exit_code resume_session=""
    local -a command=(nextflow -log "${telemetry_dir}/nextflow.log")
    if [[ -n "${database_config}" ]]; then
        command+=(-c "${database_config}")
    fi
    command+=(run "${PIPELINE_ROOT}" -profile "$(IFS=,; printf '%s' "${profiles[*]}")")
    command+=(-name "metagenomics_${run_token}_${invocation_id}" -work-dir "${invocation_work}")
    if [[ "${resume}" == true ]]; then
        resume_session="$(latest_resume_session "${invocation_id}")"
        if [[ -n "${resume_session}" ]]; then
            command+=(-resume "${resume_session}")
        else
            printf 'No prior session UUID for %s; starting this stage without -resume.\n' \
                "${invocation_id}" >&2
        fi
    fi
    command+=(--executionStage "${stage}" --outdir "${outdir}" --telemetryDir "${telemetry_dir}")
    command+=("$@")

    if [[ "${dry_run}" == true ]]; then
        print_command "${command[@]}"
        return 0
    fi

    set_storage_stage "${invocation_id}" "${stage}"
    mkdir -p -- "${telemetry_dir}" "${invocation_work}"
    started="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    set +e
    "${command[@]}"
    exit_code=$?
    set -e
    finished="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    status="completed"
    [[ ${exit_code} -eq 0 ]] || status="failed"
    record_resume_session "${invocation_id}" "${telemetry_dir}"
    append_registry "${run_token}_${invocation_id}" "${stage}" "${telemetry_dir}" \
        "${started}" "${finished}" "${status}" "${exit_code}"
    return "${exit_code}"
}

set_storage_stage() {
    local invocation_id="$1" stage="$2" temporary="${monitor_stage_file}.tmp.$$"
    [[ -n "${monitor_stage_file}" ]] || return 0
    [[ "${invocation_id}" != *$'\t'* && "${stage}" != *$'\t'* \
        && "${invocation_id}" != *$'\n'* && "${stage}" != *$'\n'* ]] \
        || die "internal storage stage context contains unsafe control characters"
    printf '%s\t%s\n' "${invocation_id}" "${stage}" > "${temporary}"
    mv -f -- "${temporary}" "${monitor_stage_file}"
}

start_storage_monitor() {
    local work_root="$1" checkpoint_root="${2:-}" cache_root="${3:-}" scratch_root="${4:-}" temp_root="${5:-}"
    local -a command=(python3 "${STORAGE_MONITOR}" --output-dir "${resource_root}/resources" \
        --work-dir "${work_root}" --results-dir "${outdir}" \
        --stop-file "${monitor_stop_file}" \
        --sample-request-file "${monitor_sample_request_file}" \
        --stage-file "${monitor_stage_file}" \
        --interval-seconds "${resource_sample_interval}")
    [[ -n "${checkpoint_root}" ]] && command+=(--checkpoint-dir "${checkpoint_root}")
    [[ -n "${cache_root}" ]] && command+=(--sra-cache-dir "${cache_root}")
    [[ -n "${scratch_root}" ]] && command+=(--sra-scratch-dir "${scratch_root}")
    [[ -n "${temp_root}" ]] && command+=(--sra-temp-dir "${temp_root}")
    [[ -n "${resource_database_root}" ]] && command+=(--database-dir "${resource_database_root}")
    set_storage_stage launcher_initializing launcher
    "${command[@]}" > "${resource_root}/resources/storage_monitor.log" 2>&1 &
    monitor_pid=$!
    telemetry_started=true
}

force_storage_sample() {
    [[ "${telemetry_started}" == true && -n "${monitor_pid}" ]] || return 0
    if ! kill -0 "${monitor_pid}" 2>/dev/null; then
        printf 'Warning: storage monitor is unavailable; cleanup will continue without a forced sample.\n' >&2
        return 0
    fi
    local token="${run_token}_$$_${RANDOM}_$(date -u +%s)" deadline=$((SECONDS + 120))
    local temporary="${monitor_sample_request_file}.tmp.$$"
    printf '%s\n' "${token}" > "${temporary}"
    mv -f -- "${temporary}" "${monitor_sample_request_file}"
    while [[ -e "${monitor_sample_request_file}" ]]; do
        if ! kill -0 "${monitor_pid}" 2>/dev/null; then
            printf 'Warning: storage monitor stopped before acknowledging the cleanup sample.\n' >&2
            return 0
        fi
        if ((SECONDS >= deadline)); then
            printf 'Warning: storage monitor did not acknowledge the cleanup sample within 120 seconds.\n' >&2
            return 0
        fi
        sleep 0.2
    done
}

summarize_resources() {
    [[ -s "${trace_registry}" ]] || return 0
    local -a command=(python3 "${RESOURCE_SUMMARIZER}" --registry "${trace_registry}" \
        --output-dir "${resource_root}" \
        --task-peaks "${resource_root}/resources/task_workdir_peaks.tsv" \
        --task-work-timeseries "${resource_root}/resources/task_workdir_timeseries.tsv" \
        --storage-timeseries "${resource_root}/resources/storage_usage_timeseries.tsv" \
        --project-status "${project_status_file}" \
        --gpu-metrics-dir "${resource_root}/resources/gpu_tasks" \
        --slurm-accounting "${resource_root}/resources/slurm_accounting.tsv")
    if [[ -n "${project_manifest}" && -f "${project_manifest}" ]]; then
        command+=(--input-manifest "${project_manifest}")
        if [[ -n "${sra_project}" ]]; then
            command+=(--manifest-mode sra)
        else
            command+=(--manifest-mode local)
        fi
    fi
    "${command[@]}" || printf 'Warning: resource summarization failed; raw telemetry was retained.\n' >&2
}

finish_telemetry() {
    local requested_exit="$1"
    trap - EXIT
    trap 'release_run_locks; exit 130' INT TERM
    set +e
    project_exit_code="${requested_exit}"
    if [[ "${project_status}" != "complete" ]]; then
        project_status="failed"
    fi
    if [[ "${telemetry_started}" == true ]]; then
        set_storage_stage launcher_finalize launcher
        : > "${monitor_stop_file}"
        wait "${monitor_pid}" 2>/dev/null || true
        rm -f -- "${monitor_stop_file}" "${monitor_sample_request_file}" "${monitor_stage_file}"
    fi
    write_project_status "${project_status}" "${project_exit_code}"
    if [[ "${environment}" == hpc ]]; then
        python3 "${SLURM_ACCOUNTING_COLLECTOR}" --registry "${trace_registry}" \
            --output "${resource_root}/resources/slurm_accounting.tsv" || true
    fi
    summarize_resources
    release_run_locks
    trap - INT TERM
    exit "${requested_exit}"
}

safe_remove_sample_work() {
    local target="$1" preprocessing_root="$2" sample_id="$3"
    local resolved_target resolved_root
    resolved_target="$(absolute_path "${target}")"
    resolved_root="$(absolute_path "${preprocessing_root}")"
    [[ "${resolved_target}" == "${resolved_root}/${sample_id}" ]] \
        || die "cleanup safety check refused sample work path: ${resolved_target}"
    [[ "${resolved_target}" != / && "${resolved_target}" != "${PIPELINE_ROOT}" ]] \
        || die "cleanup safety check refused broad path: ${resolved_target}"
    if [[ -d "${resolved_target}" ]]; then
        force_storage_sample
        rm -rf -- "${resolved_target}"
        printf 'Removed disposable completed-sample work: %s\n' "${resolved_target}"
    fi
}

safe_remove_global_work() {
    local target="${work_dir}/global" resolved_target resolved_root
    resolved_target="$(absolute_path "${target}")"
    resolved_root="$(absolute_path "${work_dir}")"
    [[ "${resolved_target}" == "${resolved_root}/global" ]] \
        || die "cleanup safety check refused global work path: ${resolved_target}"
    [[ "${resolved_target}" != / && "${resolved_target}" != "${PIPELINE_ROOT}" ]] \
        || die "cleanup safety check refused broad global work path: ${resolved_target}"
    if [[ -d "${resolved_target}" ]]; then
        force_storage_sample
        rm -rf -- "${resolved_target}"
        printf 'Removed disposable successful global work: %s\n' "${resolved_target}"
    fi
}

cleanup_sra_checkpoints() {
    local checkpoint_manifest="$1" success_marker="$2"
    local -a cleanup_arguments=(cleanup --checkpoint-manifest "${checkpoint_manifest}" \
        --checkpoint-dir "${sra_checkpoint_dir}" --success-marker "${success_marker}")
    [[ "${keep_sra_checkpoints}" == true ]] && cleanup_arguments+=(--keep)
    force_storage_sample
    python3 "${CHECKPOINT_MANAGER}" "${cleanup_arguments[@]}"
}

seal_sra_global_outputs() {
    local checkpoint_manifest="$1" success_marker="$2"
    python3 "${CHECKPOINT_MANAGER}" seal-global \
        --success-marker "${success_marker}" \
        --checkpoint-manifest "${checkpoint_manifest}" \
        --results-dir "${outdir}" \
        --project-accession "${sra_project}"
}

validate_sample_checkpoint() {
    local sample_id="$1"
    if ! python3 "${CHECKPOINT_MANAGER}" validate-sample \
        --run-manifest "${project_manifest}" \
        --checkpoint-dir "${sra_checkpoint_dir}" \
        --sample-id "${sample_id}"; then
        die "durable checkpoint validation failed for sample ${sample_id}; its work directory was retained"
    fi
}

validate_frozen_sra_state() {
    local frozen_state frozen_project frozen_platforms requested_platforms
    python3 "${SRA_RESOLVER}" --validate-existing "${state_dir}"
    frozen_state="$(python3 - "${state_dir}/sra_project_summary.json" <<'PY'
import json, sys
data = json.load(open(sys.argv[1], encoding='utf-8'))
print(data['project_accession'] + '\t' + ','.join(data['allowed_platforms']))
PY
)"
    IFS=$'\t' read -r frozen_project frozen_platforms <<< "${frozen_state}"
    requested_platforms="$(python3 - "${sra_platforms}" <<'PY'
import sys
print(','.join(sorted({value.strip().upper() for value in sys.argv[1].split(',') if value.strip()})))
PY
)"
    [[ "${frozen_project}" == "${sra_project}" ]] \
        || die "existing frozen manifest belongs to ${frozen_project}, not ${sra_project}"
    [[ "${frozen_platforms}" == "${requested_platforms}" ]] \
        || die "existing frozen manifest used platform allowlist ${frozen_platforms}; requested ${requested_platforms}; use a fresh results/state root to resolve a different cohort"
}

while (($#)); do
    case "$1" in
        --local) set_selection environment local "${environment}"; environment="local" ;;
        --hpc) set_selection environment hpc "${environment}"; environment="hpc" ;;
        --docker) set_selection runtime docker "${runtime}"; runtime="docker" ;;
        --conda) set_selection runtime conda "${runtime}"; runtime="conda" ;;
        --apptainer) set_selection runtime apptainer "${runtime}"; runtime="apptainer" ;;
        --singularity) set_selection runtime singularity "${runtime}"; runtime="singularity" ;;
        --run) set_selection mode run "${mode}"; mode="run" ;;
        --prepare-databases) set_selection mode prepare-databases "${mode}"; mode="prepare-databases" ;;
        --database-config) require_value "$1" "${2:-}"; database_config="$2"; shift ;;
        --database-config=*) database_config="${1#*=}"; [[ -n "${database_config}" ]] || die "--database-config requires a value" ;;
        --db-root) require_value "$1" "${2:-}"; db_root="$2"; shift ;;
        --db-root=*) db_root="${1#*=}" ;;
        --input) require_value "$1" "${2:-}"; input_path="$2"; forwarded_args+=(--input "$2"); shift ;;
        --input=*) input_path="${1#*=}"; [[ -n "${input_path}" ]] || die "--input requires a value"; forwarded_args+=("$1") ;;
        --sra-project) require_value "$1" "${2:-}"; sra_project="$2"; shift ;;
        --sra-project=*) sra_project="${1#*=}" ;;
        --sra-checkpoint-dir) require_value "$1" "${2:-}"; sra_checkpoint_dir="$2"; shift ;;
        --sra-checkpoint-dir=*) sra_checkpoint_dir="${1#*=}" ;;
        --sra-scratch-dir) require_value "$1" "${2:-}"; sra_scratch_root="$2"; shift ;;
        --sra-scratch-dir=*) sra_scratch_root="${1#*=}" ;;
        --sra-cache-dir) require_value "$1" "${2:-}"; sra_cache_dir="$2"; shift ;;
        --sra-cache-dir=*) sra_cache_dir="${1#*=}" ;;
        --sra-temp-dir) require_value "$1" "${2:-}"; sra_temp_dir="$2"; shift ;;
        --sra-temp-dir=*) sra_temp_dir="${1#*=}" ;;
        --sra-email) require_value "$1" "${2:-}"; sra_email="$2"; shift ;;
        --sra-email=*) sra_email="${1#*=}" ;;
        --sra-platforms) require_value "$1" "${2:-}"; sra_platforms="$2"; shift ;;
        --sra-platforms=*) sra_platforms="${1#*=}" ;;
        --sra-max-size) require_value "$1" "${2:-}"; sra_max_size="$2"; shift ;;
        --sra-max-size=*) sra_max_size="${1#*=}" ;;
        --keep-sra-checkpoints) keep_sra_checkpoints=true ;;
        --outdir) require_value "$1" "${2:-}"; outdir="$2"; shift ;;
        --outdir=*) outdir="${1#*=}" ;;
        --work-dir) require_value "$1" "${2:-}"; work_dir="$2"; shift ;;
        --work-dir=*) work_dir="${1#*=}" ;;
        --publish_dir_mode) require_value "$1" "${2:-}"; publish_dir_mode_override="$2"; forwarded_args+=("$1" "$2"); shift ;;
        --publish_dir_mode=*) publish_dir_mode_override="${1#*=}"; [[ -n "${publish_dir_mode_override}" ]] || die "--publish_dir_mode requires a value"; forwarded_args+=("$1") ;;
        --resource-database-root) require_value "$1" "${2:-}"; resource_database_root="$2"; shift ;;
        --resource-database-root=*) resource_database_root="${1#*=}" ;;
        --resource-sample-interval) require_value "$1" "${2:-}"; resource_sample_interval="$2"; shift ;;
        --resource-sample-interval=*) resource_sample_interval="${1#*=}" ;;
        --storage-constrained) storage_constrained=true ;;
        --enable-gpu) enable_gpu=true ;;
        --gpu-accelerators) require_value "$1" "${2:-}"; gpu_accelerators="$2"; shift ;;
        --gpu-accelerators=*) gpu_accelerators="${1#*=}" ;;
        --gpu-telemetry-interval) require_value "$1" "${2:-}"; gpu_telemetry_interval="$2"; shift ;;
        --gpu-telemetry-interval=*) gpu_telemetry_interval="${1#*=}" ;;
        --slurm-gpu-gres) require_value "$1" "${2:-}"; slurm_gpu_gres="$2"; shift ;;
        --slurm-gpu-gres=*) slurm_gpu_gres="${1#*=}" ;;
        --resume) resume=true ;;
        --dry-run) dry_run=true ;;
        --version)
            printf 'metagenomics_pipeline.sh 2.0.0\nRequired Nextflow: >=%s\n' "${MINIMUM_NEXTFLOW_VERSION}"
            exit 0
            ;;
        -h|--help) usage; exit 0 ;;
        --) die "the '--' passthrough delimiter is not supported; pass pipeline parameters directly" ;;
        -profile|-profile=*|-work-dir|-work-dir=*|-w|-w=*|-name|-name=*|-log|-log=*|-resume|-resume=*|-params-file|-params-file=*|-c|-c=*|-config|-config=*)
            die "${1%%=*} is reserved by the staged launcher"
            ;;
        --executionStage|--executionStage=*|--telemetryDir|--telemetryDir=*|--sraStateDir|--sraStateDir=*|--sraManifest|--sraManifest=*|--sraSampleId|--sraSampleId=*|--sraCheckpointManifest|--sraCheckpointManifest=*|--sraRequireComplete|--sraRequireComplete=*|--sraContainerOptions|--sraContainerOptions=*|--sraProject|--sraProject=*|--sraCheckpointDir|--sraCheckpointDir=*|--sraScratchDir|--sraScratchDir=*|--sraCacheDir|--sraCacheDir=*|--sraTempDir|--sraTempDir=*|--enableGpu|--enableGpu=*|--gpuAccelerators|--gpuAccelerators=*|--gpuTelemetryInterval|--gpuTelemetryInterval=*|--gpuContainerOptions|--gpuContainerOptions=*|--slurmGpuGres|--slurmGpuGres=*)
            die "${1%%=*} is an internal pipeline parameter managed by the launcher"
            ;;
        *) forwarded_args+=("$1") ;;
    esac
    shift
done

[[ -n "${environment}" ]] || die "select --local or --hpc"
[[ -n "${mode}" ]] || die "select one execution mode"

if [[ "${mode}" == "prepare-databases" ]]; then
    [[ -z "${runtime}" ]] || die "software runtime flags do not apply to --prepare-databases"
    [[ -n "${db_root}" ]] || die "--prepare-databases requires --db-root PATH"
    [[ -f "${DATABASE_PREPARER}" ]] || die "database preparation script was not found"
    declare -a prepare_command=(bash "${DATABASE_PREPARER}" --db-root "${db_root}")
    [[ "${dry_run}" == true ]] && prepare_command+=(--dry-run)
    prepare_command+=("${forwarded_args[@]}")
    if [[ "${dry_run}" == true ]]; then
        print_command "${prepare_command[@]}"
    else
        "${prepare_command[@]}"
    fi
    exit 0
fi

[[ -n "${runtime}" ]] || die "select one software runtime"
[[ -z "${db_root}" ]] || die "--db-root only applies to --prepare-databases"
[[ -n "${input_path}" || -n "${sra_project}" ]] || die "--run requires exactly one of --input or --sra-project"
[[ -z "${input_path}" || -z "${sra_project}" ]] || die "--input and --sra-project are mutually exclusive"
[[ "${environment}" != hpc || "${runtime}" != docker ]] || die "Docker is not supported by the SLURM launcher"
[[ -z "${database_config}" || -r "${database_config}" ]] || die "database configuration is not readable: ${database_config}"
[[ "${gpu_accelerators}" =~ ^[1-9][0-9]*$ ]] || die "--gpu-accelerators must be a positive integer"
[[ "${gpu_accelerators}" == 1 ]] \
    || die "the verified COMEBin/SemiBin2/Vamb paths each use exactly one GPU; --gpu-accelerators must be 1"
[[ "${gpu_telemetry_interval}" =~ ^[1-9][0-9]*$ ]] \
    || die "--gpu-telemetry-interval must be a positive integer"
[[ "${resource_sample_interval}" =~ ^(0*[1-9][0-9]*|[0-9]*\.[0-9]*[1-9][0-9]*)$ ]] \
    || die "--resource-sample-interval must be greater than zero"
[[ -z "${slurm_gpu_gres}" || "${slurm_gpu_gres}" =~ ^[A-Za-z0-9_.:+,-]+$ ]] \
    || die "--slurm-gpu-gres contains unsupported characters"
[[ "${environment}" != hpc || "${enable_gpu}" != true || -n "${slurm_gpu_gres}" ]] \
    || die "HPC GPU mode requires an explicit --slurm-gpu-gres value"
if [[ "${environment}" == local && "${enable_gpu}" == true && "${runtime}" != docker ]]; then
    local_cuda_device="${CUDA_VISIBLE_DEVICES:-}"
    [[ -n "${local_cuda_device}" && "${local_cuda_device}" != *,* \
        && "${local_cuda_device}" != *[[:space:]]* ]] \
        || die "local Conda/Apptainer/Singularity GPU mode requires CUDA_VISIBLE_DEVICES to name exactly one device"
fi

profiles=("${environment}" "${runtime}")
if [[ "${storage_constrained}" == true || -n "${sra_project}" || ( "${enable_gpu}" == true && "${environment}" == local ) ]]; then
    # The local executor does not schedule the accelerator directive.  A queue
    # size of one prevents independently ready GPU tools from contending for
    # the same device; SRA also serializes independent global branches to keep
    # their large temporary trees from overlapping.  SLURM GPU allocation uses GRES.
    profiles+=(disk_efficient)
fi
[[ "${enable_gpu}" == true ]] && profiles+=(gpu)

if [[ "${dry_run}" == false ]]; then
    check_nextflow
    check_runtime
fi

outdir="$(absolute_path "${outdir}")"
if [[ -n "${resource_database_root}" ]]; then
    resource_database_root="$(absolute_path "${resource_database_root}")"
fi
resource_root="${outdir}/pipeline_info"
trace_registry="${resource_root}/resources/trace_registry.tsv"
resume_registry="${resource_root}/resources/resume_sessions.tsv"
project_status_file="${resource_root}/resources/project_status.json"
run_token="$(date -u +%Y%m%dT%H%M%SZ)_$$"
monitor_stop_file="${resource_root}/resources/.storage_monitor_stop_${run_token}"
monitor_sample_request_file="${resource_root}/resources/.storage_sample_request_${run_token}"
monitor_stage_file="${resource_root}/resources/.storage_stage_${run_token}.tsv"

if [[ "${dry_run}" == false ]]; then
    check_python
    mkdir -p -- "${resource_root}"
    run_lock_token="${run_token}_${RANDOM}_${RANDOM}"
    trap 'release_run_locks' EXIT
    trap 'exit 130' INT TERM
    acquire_run_lock "${resource_root}/.metagenomics_run.lock" results
fi

declare -a common_nextflow_args=("${forwarded_args[@]}")
common_nextflow_args+=(--enableGpu "${enable_gpu}" --gpuAccelerators "${gpu_accelerators}" \
    --gpuTelemetryInterval "${gpu_telemetry_interval}")
[[ -n "${slurm_gpu_gres}" ]] && common_nextflow_args+=(--slurmGpuGres "${slurm_gpu_gres}")
if [[ "${enable_gpu}" == true && "${runtime}" == docker ]]; then
    common_nextflow_args+=(--gpuContainerOptions '--gpus 1')
elif [[ "${enable_gpu}" == true && ( "${runtime}" == apptainer || "${runtime}" == singularity ) ]]; then
    common_nextflow_args+=(--gpuContainerOptions '--nv')
fi

if [[ -n "${input_path}" ]]; then
    project_manifest="$(absolute_path "${input_path}")"
    work_dir="$(absolute_path "${work_dir:-${PWD}/work}")"
    if [[ "${dry_run}" == false ]]; then
        initialize_project_telemetry
        mkdir -p -- "${work_dir}" "${outdir}"
        start_storage_monitor "${work_dir}" "" "" "" ""
    fi
    run_nextflow local local "${work_dir}" "${common_nextflow_args[@]}"
    project_status="complete"
    exit 0
fi

[[ "${sra_project}" =~ ^PRJ(NA|EB|DB)[0-9]+$ ]] || die "invalid BioProject accession: ${sra_project}"
[[ -z "${publish_dir_mode_override}" || "${publish_dir_mode_override}" == copy ]] \
    || die "SRA mode requires --publish_dir_mode copy because transient work is cleaned after success"
[[ -n "${sra_checkpoint_dir}" ]] || die "SRA mode requires --sra-checkpoint-dir outside transient work"
[[ -n "${sra_scratch_root}" ]] || die "SRA mode requires --sra-scratch-dir"
[[ "${sra_platforms}" =~ ^[A-Za-z0-9_.-]+(,[A-Za-z0-9_.-]+)*$ ]] \
    || die "--sra-platforms must be a comma-separated list of platform names"
[[ "${sra_max_size}" =~ ^(u|[1-9][0-9]*([KMGT]B?)?)$ ]] \
    || die "--sra-max-size must be u or a positive size accepted by prefetch"
[[ -z "${sra_email}" || "${sra_email}" =~ ^[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+$ ]] \
    || die "--sra-email contains unsupported characters or is not an email address"
sra_checkpoint_dir="$(absolute_path "${sra_checkpoint_dir}")"
sra_scratch_root="$(absolute_path "${sra_scratch_root}")"
sra_cache_dir="$(absolute_path "${sra_cache_dir:-${sra_scratch_root}/sra-cache}")"
sra_temp_dir="$(absolute_path "${sra_temp_dir:-${sra_scratch_root}/fasterq-temp}")"
work_dir="$(absolute_path "${work_dir:-${sra_scratch_root}/nextflow-work}")"
acquisition_scratch="${sra_scratch_root}/acquisition"
sra_container_home="${acquisition_scratch}/container-home"
preprocessing_work_root="${work_dir}/preprocess"
state_dir="${resource_root}/sra"
project_manifest="${state_dir}/sra_project_manifest.tsv"
checkpoint_manifest="${state_dir}/sra_checkpoint_manifest.tsv"

for path in "${sra_checkpoint_dir}" "${sra_scratch_root}" "${sra_cache_dir}" "${sra_temp_dir}" "${work_dir}"; do
    [[ "${path}" != / && "${path}" != "${PIPELINE_ROOT}" ]] || die "unsafe SRA storage path: ${path}"
    require_interpolation_safe_path "SRA paths" "${path}"
done
[[ "${outdir}" != / ]] || die "unsafe results path: ${outdir}"
require_interpolation_safe_path "--outdir" "${outdir}"

declare -a disposable_roots=("${sra_scratch_root}" "${sra_cache_dir}" "${sra_temp_dir}" "${work_dir}")
for path in "${sra_checkpoint_dir}" "${disposable_roots[@]}"; do
    path_is_within "${path}" "${PIPELINE_ROOT}" \
        && die "SRA checkpoints, cache, temporary data, scratch, and work must be outside the Git repository: ${path}"
done
for path in "${disposable_roots[@]}"; do
    if path_is_within "${sra_checkpoint_dir}" "${path}" \
        || path_is_within "${path}" "${sra_checkpoint_dir}"; then
        die "the durable SRA checkpoint root must not overlap disposable storage: ${path}"
    fi
    if path_is_within "${outdir}" "${path}" \
        || path_is_within "${path}" "${outdir}"; then
        die "the results root must not overlap disposable SRA storage: ${path}"
    fi
done
if path_is_within "${outdir}" "${sra_checkpoint_dir}" \
    || path_is_within "${sra_checkpoint_dir}" "${outdir}"; then
    die "the results and durable SRA checkpoint roots must not overlap"
fi
declare -a measured_disposable_roots=(
    "${acquisition_scratch}" "${sra_cache_dir}" "${sra_temp_dir}" "${work_dir}"
)
for ((left_index = 0; left_index < ${#measured_disposable_roots[@]}; left_index++)); do
    for ((right_index = left_index + 1; right_index < ${#measured_disposable_roots[@]}; right_index++)); do
        left_path="${measured_disposable_roots[left_index]}"
        right_path="${measured_disposable_roots[right_index]}"
        if path_is_within "${left_path}" "${right_path}" \
            || path_is_within "${right_path}" "${left_path}"; then
            die "separately measured SRA cache, temp, acquisition, and work roots must not overlap: ${left_path} / ${right_path}"
        fi
    done
done
if [[ "${runtime}" == docker || "${runtime}" == apptainer || "${runtime}" == singularity ]]; then
    for path in "${sra_checkpoint_dir}" "${sra_scratch_root}" "${sra_cache_dir}" "${sra_temp_dir}"; do
        require_container_safe_path "SRA storage paths" "${path}"
    done
fi

if [[ "${dry_run}" == false ]]; then
    mkdir -p -- "$(dirname -- "${sra_checkpoint_dir}")"
    acquire_run_lock "${sra_checkpoint_dir}.metagenomics_run.lock" checkpoint \
        "${sra_checkpoint_dir}"
    initialize_project_telemetry
fi

if [[ "${runtime}" == docker ]]; then
    sra_container_options="--env HOME=${sra_container_home} -v ${sra_checkpoint_dir}:${sra_checkpoint_dir} -v ${sra_scratch_root}:${sra_scratch_root} -v ${sra_cache_dir}:${sra_cache_dir} -v ${sra_temp_dir}:${sra_temp_dir}"
elif [[ "${runtime}" == apptainer || "${runtime}" == singularity ]]; then
    sra_container_options="--bind ${sra_checkpoint_dir}:${sra_checkpoint_dir},${sra_scratch_root}:${sra_scratch_root},${sra_cache_dir}:${sra_cache_dir},${sra_temp_dir}:${sra_temp_dir}"
else
    sra_container_options=""
fi

declare -a sra_common=("${common_nextflow_args[@]}" \
    --sraProject "${sra_project}" --sraCheckpointDir "${sra_checkpoint_dir}" \
    --sraScratchDir "${acquisition_scratch}" --sraCacheDir "${sra_cache_dir}" \
    --sraTempDir "${sra_temp_dir}" --sraStateDir "${state_dir}" \
    --sraPlatforms "${sra_platforms}" --sraMaxSize "${sra_max_size}" \
    --sraContainerOptions "${sra_container_options}" \
    --save_clean_reads false --save_host_removed_reads false \
    --publish_dir_mode copy)
[[ -n "${sra_email}" ]] && sra_common+=(--sraEmail "${sra_email}")

if [[ "${dry_run}" == true ]]; then
    run_nextflow discovery sra-discovery "${work_dir}/discovery" "${sra_common[@]}"
    run_nextflow checkpoints_initial sra-checkpoints "${work_dir}/checkpoints-initial" \
        "${sra_common[@]}" --sraManifest "${project_manifest}" --sraRequireComplete false
    run_nextflow sample_SAMPLE_ID sra-preprocess "${preprocessing_work_root}/SAMPLE_ID" \
        "${sra_common[@]}" --sraManifest "${project_manifest}" --sraSampleId SAMPLE_ID
    run_nextflow checkpoints_final sra-checkpoints "${work_dir}/checkpoints-final" \
        "${sra_common[@]}" --sraManifest "${project_manifest}" --sraRequireComplete true
    run_nextflow global sra-global "${work_dir}/global" "${sra_common[@]}" \
        --sraCheckpointManifest "${state_dir}/sra_checkpoint_manifest.tsv"
    printf 'Dry run only: SAMPLE_ID represents the deterministic pending-sample loop.\n'
    exit 0
fi

mkdir -p -- "${sra_checkpoint_dir}" "${acquisition_scratch}" "${sra_cache_dir}" \
    "${sra_temp_dir}" "${sra_container_home}" "${preprocessing_work_root}" "${state_dir}" "${outdir}"
start_storage_monitor "${work_dir}" "${sra_checkpoint_dir}" "${sra_cache_dir}" \
    "${acquisition_scratch}" "${sra_temp_dir}"

if [[ -f "${project_manifest}" ]]; then
    validate_frozen_sra_state
fi

if [[ -f "${state_dir}/sra_global_success.json" ]]; then
    [[ -f "${project_manifest}" ]] \
        || die "global-success marker exists but the frozen SRA manifest is missing"
    [[ -f "${checkpoint_manifest}" ]] \
        || die "global-success marker exists but the checkpoint manifest is missing"
    if ! python3 - "${state_dir}/sra_global_success.json" "${sra_project}" <<'PY'
import json, sys
data = json.load(open(sys.argv[1], encoding='utf-8'))
raise SystemExit(0 if data.get('status') == 'complete' and data.get('project_accession') == sys.argv[2] else 1)
PY
    then
        die "global-success marker does not match ${sra_project} or is not complete"
    fi
    seal_sra_global_outputs "${checkpoint_manifest}" "${state_dir}/sra_global_success.json"
    cleanup_sra_checkpoints "${checkpoint_manifest}" "${state_dir}/sra_global_success.json"
    safe_remove_global_work
    printf 'Project %s already has validated durable global outputs; no scientific stage was rerun.\n' "${sra_project}"
    project_status="complete"
    exit 0
fi

if [[ ! -f "${project_manifest}" ]]; then
    run_nextflow discovery sra-discovery "${work_dir}/discovery" "${sra_common[@]}"
    validate_frozen_sra_state
fi

run_nextflow checkpoints_initial sra-checkpoints "${work_dir}/checkpoints-initial" \
    "${sra_common[@]}" --sraManifest "${project_manifest}" --sraRequireComplete false

mapfile -t pending_samples < <(awk -F '\t' 'NR > 1 && $2 != "" { print $2 }' "${state_dir}/sra_pending_samples.tsv")
for sample_id in "${pending_samples[@]}"; do
    [[ "${sample_id}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || die "unsafe sample ID in pending manifest: ${sample_id}"
    printf 'Processing SRA biological sample %s\n' "${sample_id}"
    sample_work="${preprocessing_work_root}/${sample_id}"
    run_nextflow "sample_${sample_id}" sra-preprocess "${sample_work}" \
        "${sra_common[@]}" --sraManifest "${project_manifest}" --sraSampleId "${sample_id}"
    validate_sample_checkpoint "${sample_id}"
    safe_remove_sample_work "${sample_work}" "${preprocessing_work_root}" "${sample_id}"
done

run_nextflow checkpoints_final sra-checkpoints "${work_dir}/checkpoints-final" \
    "${sra_common[@]}" --sraManifest "${project_manifest}" --sraRequireComplete true

run_nextflow global sra-global "${work_dir}/global" "${sra_common[@]}" \
    --sraCheckpointManifest "${checkpoint_manifest}"

success_marker="${state_dir}/sra_global_success.json"
[[ -s "${success_marker}" ]] || die "global workflow exited without its validated success marker"
seal_sra_global_outputs "${checkpoint_manifest}" "${success_marker}"
cleanup_sra_checkpoints "${checkpoint_manifest}" "${success_marker}"
safe_remove_global_work

project_status="complete"
exit 0
