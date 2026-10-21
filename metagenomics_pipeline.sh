#!/usr/bin/env bash
# Argument adapter for ONE Nextflow execution; Nextflow owns sample scheduling.
set -euo pipefail
readonly PIPELINE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
environment="" runtime="" mode="" database_config="" db_root=""
outdir="results" work_dir="${PWD}/work" dry_run=false resume=false
storage_constrained=false enable_gpu=false
input_path="" sra_project="" sra_samples="" slurm_gpu_gres=""
lock_dir="" run_token="$(date -u +%Y%m%dT%H%M%SZ)_$$"
declare -a args=() profiles=()

usage() {
    cat <<'USAGE'
Usage:
  ./metagenomics_pipeline.sh --hpc --apptainer --run --sra-project PRJNA... --sra-samples FILE [options]
  ./metagenomics_pipeline.sh --local --docker --run --input FILE [options]
  ./metagenomics_pipeline.sh --local --prepare-databases --db-root PATH [options]

Select one environment: --local | --hpc
Select one runtime for --run: --docker | --conda | --apptainer | --singularity
Input: --input CSV OR --sra-project ACCESSION plus --sra-samples TSV
Options:
  --database-config FILE       Existing database/site Nextflow configuration.
  --outdir PATH                Published results (default: results).
  --work-dir PATH              Persistent Nextflow work directory.
  --group-column NAME          Group metadata column; group is detected by default.
  --spades-coassembly-mode MODE  global (default) or condition.
  --sra-run-info CSV           Existing NCBI RunInfo metadata (optional).
  --sra-download-max-forks N   Maximum simultaneous SRA run tasks (default: 4).
  --sra-scratch-dir PATH       Optional Nextflow scratch root for SRA tasks.
  --sra-platforms LIST         Allowed platforms (default: ILLUMINA,BGISEQ).
  --sra-max-size VALUE         Per-run prefetch limit (default: u).
  --sra-email EMAIL            NCBI metadata contact email.
  --storage-constrained       Limit acquisition concurrency and queued tasks;
                              keep cached outputs so -resume remains usable.
  --enable-gpu                Existing COMEBin/SemiBin2/Vamb GPU profile.
  --gpu-accelerators N        One GPU per task (only 1 is supported).
  --gpu-telemetry-interval N  GPU sampling interval (default: 10 seconds).
  --slurm-gpu-gres VALUE      Explicit site-specific SLURM GRES.
  --resume                    Forward native Nextflow -resume (same launch directory).
  --dry-run                   Print the single command without writes or execution.
  --version | --help

Other pipeline parameters are forwarded to Nextflow. No sample loops, external
completion records, work cleanup, or separate workflow invocations are used.
USAGE
}
die() { printf 'Error: %s\n' "$*" >&2; exit 2; }
require_value() { [[ -n "${2:-}" && "$2" != --* ]] || die "$1 requires a value"; }
select_one() {
    local previous="$1" selected="$2" kind="$3"
    [[ -z "$previous" || "$previous" == "$selected" ]] || die "select exactly one $kind"
}
absolute_path() { realpath -m -- "$1"; }
print_command() { printf 'Command:'; printf ' %q' "$@"; printf '\n'; }

while (($#)); do
    if [[ "$1" == --*=* ]]; then
        option="${1%%=*}"; value="${1#*=}"; shift
        case "$option" in
            --local|--hpc|--docker|--conda|--apptainer|--singularity|--run|--prepare-databases|--storage-constrained|--enable-gpu|--resume|--dry-run)
                die "$option is a flag and takes no value" ;;
        esac
        set -- "$option" "$value" "$@"
    fi
    case "$1" in
        --local|--hpc) select_one "$environment" "${1#--}" environment; environment="${1#--}" ;;
        --docker|--conda|--apptainer|--singularity) select_one "$runtime" "${1#--}" runtime; runtime="${1#--}" ;;
        --run|--prepare-databases) select_one "$mode" "${1#--}" mode; mode="${1#--}" ;;
        --database-config|--outdir|--work-dir|--db-root)
            option="$1"; require_value "$1" "${2:-}"; value="$2"; shift
            case "$option" in
                --database-config) database_config="$value" ;;
                --outdir) outdir="$value" ;;
                --work-dir) work_dir="$value" ;;
                --db-root) db_root="$value" ;;
            esac ;;
        --input|--sra-project|--sra-samples|--group-column|--spades-coassembly-mode|--sra-run-info|--sra-download-max-forks|--sra-scratch-dir|--sra-platforms|--sra-max-size|--sra-email|--gpu-accelerators|--gpu-telemetry-interval|--slurm-gpu-gres)
            option="$1"; require_value "$1" "${2:-}"; value="$2"; shift
            case "$option" in
                --input) input_path="$value"; key=input ;;
                --sra-project) sra_project="$value"; key=sraProject ;;
                --sra-samples) sra_samples="$value"; key=sraSamples ;;
                --group-column) key=groupColumn ;;
                --spades-coassembly-mode) key=spadesCoassemblyMode ;;
                --sra-run-info) key=sraRunInfo ;;
                --sra-download-max-forks) key=sraDownloadMaxForks ;;
                --sra-scratch-dir) key=sraScratchDir ;;
                --sra-platforms) key=sraPlatforms ;;
                --sra-max-size) key=sraMaxSize ;;
                --sra-email) key=sraEmail ;;
                --gpu-accelerators) [[ "$value" == 1 ]] || die 'only one GPU per task is supported'; key=gpuAccelerators ;;
                --gpu-telemetry-interval) [[ "$value" =~ ^[1-9][0-9]*$ ]] || die 'GPU interval must be positive'; key=gpuTelemetryInterval ;;
                --slurm-gpu-gres) slurm_gpu_gres="$value"; key=slurmGpuGres ;;
            esac
            args+=("--$key" "$value") ;;
        --storage-constrained) storage_constrained=true ;;
        --enable-gpu) enable_gpu=true ;;
        --resume) resume=true ;;
        --dry-run) dry_run=true ;;
        --version) printf 'metagenomics_pipeline.sh 3.0.0\nNextflow >=26.04.6\n'; exit 0 ;;
        -h|--help) usage; exit 0 ;;
        -resume|-profile|-work-dir|-w|-c|-config|-log|-name|--)
            die "$1: use the corresponding launcher option or invoke Nextflow directly" ;;
        *) args+=("$1") ;;
    esac
    shift
done
[[ -n "$environment" && -n "$mode" ]] || die 'select an environment and an execution mode'
if [[ "$mode" == prepare-databases ]]; then
    [[ -z "$runtime" && -n "$db_root" ]] || die 'database preparation requires --db-root and no runtime flag'
    command=(bash "$PIPELINE_ROOT/bin/prepare_databases.sh" --db-root "$db_root" "${args[@]}")
    if [[ "$dry_run" == true ]]; then print_command "${command[@]}"; else exec "${command[@]}"; fi
    exit 0
fi
[[ -n "$runtime" && -z "$db_root" ]] || die 'select a runtime; --db-root is only for database preparation'
[[ -n "$input_path" || ( -n "$sra_project" && -n "$sra_samples" ) ]] || die 'provide --input or --sra-project plus --sra-samples'
[[ -z "$input_path" || ( -z "$sra_project" && -z "$sra_samples" ) ]] || die 'input modes are mutually exclusive'
[[ -z "$database_config" || -r "$database_config" ]] || die 'database configuration is not readable'
[[ "$environment" != hpc || "$runtime" != docker ]] || die 'Docker is not supported by this SLURM launcher'
profiles=("$environment" "$runtime")
[[ "$storage_constrained" == true ]] && profiles+=(disk_efficient)
if [[ "$enable_gpu" == true ]]; then
    profiles+=(gpu)
    args+=(--enableGpu true)
    if [[ "$environment" == hpc ]]; then
        [[ -n "$slurm_gpu_gres" && "$slurm_gpu_gres" =~ ^[A-Za-z0-9_.:+,-]+$ ]] || die 'HPC GPU mode requires --slurm-gpu-gres'
    else
        # The local executor has no GPU allocator. Preserve exclusive device use.
        [[ "$storage_constrained" == true ]] || profiles+=(disk_efficient)
        args+=(--storage_max_jobs 1)
        if [[ "$runtime" != docker ]]; then
            [[ -n "${CUDA_VISIBLE_DEVICES:-}" && "$CUDA_VISIBLE_DEVICES" != *,* && "$CUDA_VISIBLE_DEVICES" != *[[:space:]]* ]] || die 'select exactly one CUDA_VISIBLE_DEVICES device'
        fi
    fi
    case "$runtime" in
        docker) args+=(--gpuContainerOptions '--gpus 1') ;;
        apptainer|singularity) args+=(--gpuContainerOptions '--nv') ;;
    esac
fi
outdir="$(absolute_path "$outdir")"
work_dir="$(absolute_path "$work_dir")"
info="$outdir/pipeline_info"
invocation="$info/invocations/$run_token"
command=(nextflow -log "$invocation/nextflow.log")
[[ -z "$database_config" ]] || command+=(-c "$(absolute_path "$database_config")")
command+=(run "$PIPELINE_ROOT" -profile "$(IFS=,; printf '%s' "${profiles[*]}")" -work-dir "$work_dir")
[[ "$resume" == false ]] || command+=(-resume)
command+=(--outdir "$outdir" --reportsDir "$invocation" "${args[@]}")
if [[ "$dry_run" == true ]]; then print_command "${command[@]}"; exit 0; fi
command -v nextflow >/dev/null || die 'Nextflow is not available on PATH'
case "$runtime" in
    conda) command -v mamba >/dev/null || die 'Mamba is required' ;;
    docker) docker info >/dev/null 2>&1 || die 'Docker daemon is unavailable' ;;
    *) command -v "$runtime" >/dev/null || die "$runtime is not available" ;;
esac
mkdir -p -- "$info"
lock_dir="$info/.metagenomics_run.lock"
mkdir -- "$lock_dir" 2>/dev/null || die "results lock exists: $lock_dir; inspect it before removing a stale lock"
printf '%s\n' "$run_token" > "$lock_dir/owner"
finish() {
    local code="$?"
    trap - EXIT
    set +e
    if [[ -f "$lock_dir/owner" && "$(cat "$lock_dir/owner")" == "$run_token" ]]; then
        rm -- "$lock_dir/owner"
        rmdir -- "$lock_dir"
    fi
    exit "$code"
}
trap finish EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
mkdir -p -- "$invocation" "$work_dir"
"${command[@]}"
