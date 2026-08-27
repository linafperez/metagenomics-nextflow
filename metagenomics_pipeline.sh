#!/usr/bin/env bash

set -euo pipefail

readonly MINIMUM_NEXTFLOW_VERSION="26.04.6"
readonly PIPELINE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
readonly DATABASE_PREPARER="${PIPELINE_ROOT}/bin/prepare_databases.sh"

environment=""
runtime=""
mode=""
database_config=""
db_root=""
dry_run=false
resume=false
declare -a forwarded_args=()

usage() {
    cat <<'USAGE'
Usage:
  ./metagenomics_pipeline.sh --<environment> --<runtime> --<mode> [options]
  ./metagenomics_pipeline.sh --<environment> --prepare-databases --db-root PATH [options]

Environment (select one):
  --local                  Use the local executor.
  --hpc                    Use the SLURM executor.

Software runtime (select one for pipeline execution):
  --docker                 Use Docker containers.
  --conda                  Use Conda environments.
  --apptainer              Use Apptainer containers.
  --singularity            Use the Singularity compatibility profile.

Mode (select one):
  --run                    Run the complete production workflow.
  --prepare-databases      Prepare shared production databases separately.

Launcher options:
  --database-config FILE   Add a generated database configuration with Nextflow -c.
  --db-root PATH           Database root for --prepare-databases.
  --resume                 Resume the Nextflow run.
  --dry-run                Print commands without executing them.
  -h, --help               Show this help.
  --version                Show the launcher and required Nextflow versions.

Arguments not recognized by the launcher are passed unchanged to Nextflow or,
in database-preparation mode, to prepare_databases.sh. Use -- to end launcher
option processing explicitly.

Examples:
  ./metagenomics_pipeline.sh --local --docker --run --input samplesheet.csv
  ./metagenomics_pipeline.sh --hpc --apptainer --run \
      --database-config /shared/db/metagenomics_databases.config \
      --input samplesheet.hpc.csv
  ./metagenomics_pipeline.sh --hpc --prepare-databases \
      --db-root /shared/databases/metagenomics
USAGE
}

die() {
    printf 'Error: %s\n' "$*" >&2
    exit 2
}

set_selection() {
    local selection_name="$1"
    local selection_value="$2"
    local current_value="$3"

    if [[ -n "${current_value}" && "${current_value}" != "${selection_value}" ]]; then
        die "select exactly one ${selection_name}"
    fi
}

require_value() {
    local option="$1"
    local value="${2:-}"
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

version_at_least() {
    local actual="$1"
    local required="$2"
    local actual_part required_part index
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
        if ((10#${actual_part} > 10#${required_part})); then
            return 0
        fi
        if ((10#${actual_part} < 10#${required_part})); then
            return 1
        fi
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

check_runtime() {
    case "${runtime}" in
        docker)
            command -v docker >/dev/null 2>&1 || die "Docker is not available on PATH"
            docker info >/dev/null 2>&1 || die "the Docker daemon is not available"
            ;;
        conda)
            command -v mamba >/dev/null 2>&1 \
                || die "Mamba is required by the Conda profile and is not available on PATH"
            ;;
        apptainer)
            command -v apptainer >/dev/null 2>&1 || die "Apptainer is not available on PATH"
            ;;
        singularity)
            command -v singularity >/dev/null 2>&1 || die "Singularity is not available on PATH"
            ;;
        *)
            die "unsupported software runtime '${runtime}'"
            ;;
    esac
}

validate_forwarded_input() {
    local index argument next_argument
    for ((index = 0; index < ${#forwarded_args[@]}; index++)); do
        argument="${forwarded_args[index]}"
        case "${argument}" in
            --input=*)
                [[ -n "${argument#*=}" ]] || die "--input requires a samplesheet path"
                return 0
                ;;
            --input)
                next_argument="${forwarded_args[index + 1]:-}"
                [[ -n "${next_argument}" && "${next_argument}" != --* ]] \
                    || die "--input requires a samplesheet path"
                return 0
                ;;
        esac
    done
    die "${mode} requires --input PATH"
}

while (($#)); do
    case "$1" in
        --local)
            set_selection environment local "${environment}"
            environment="local"
            ;;
        --hpc)
            set_selection environment hpc "${environment}"
            environment="hpc"
            ;;
        --docker)
            set_selection runtime docker "${runtime}"
            runtime="docker"
            ;;
        --conda)
            set_selection runtime conda "${runtime}"
            runtime="conda"
            ;;
        --apptainer)
            set_selection runtime apptainer "${runtime}"
            runtime="apptainer"
            ;;
        --singularity)
            set_selection runtime singularity "${runtime}"
            runtime="singularity"
            ;;
        --run)
            set_selection mode run "${mode}"
            mode="run"
            ;;
        --prepare-databases)
            set_selection mode prepare-databases "${mode}"
            mode="prepare-databases"
            ;;
        --database-config)
            require_value "$1" "${2:-}"
            database_config="$2"
            shift
            ;;
        --database-config=*)
            database_config="${1#*=}"
            [[ -n "${database_config}" ]] || die "--database-config requires a value"
            ;;
        --db-root)
            require_value "$1" "${2:-}"
            db_root="$2"
            shift
            ;;
        --db-root=*)
            db_root="${1#*=}"
            [[ -n "${db_root}" ]] || die "--db-root requires a value"
            ;;
        --resume)
            resume=true
            ;;
        --dry-run)
            dry_run=true
            ;;
        --version)
            printf 'metagenomics_pipeline.sh 1.0.0\n'
            printf 'Required Nextflow: >=%s\n' "${MINIMUM_NEXTFLOW_VERSION}"
            exit 0
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        --)
            shift
            forwarded_args+=("$@")
            break
            ;;
        *)
            forwarded_args+=("$1")
            ;;
    esac
    shift
done

[[ -n "${environment}" ]] || die "select --local or --hpc"
[[ -n "${mode}" ]] || die "select one execution mode"

if [[ "${mode}" == "prepare-databases" ]]; then
    [[ -z "${runtime}" ]] \
        || die "software runtime flags do not apply to --prepare-databases"
    [[ -z "${database_config}" ]] \
        || die "--database-config does not apply to --prepare-databases"
    [[ "${resume}" == false ]] \
        || die "--resume does not apply to --prepare-databases"
    [[ -n "${db_root}" ]] || die "--prepare-databases requires --db-root PATH"
    [[ -f "${DATABASE_PREPARER}" ]] \
        || die "database preparation script was not found: ${DATABASE_PREPARER}"

    prepare_command=(bash "${DATABASE_PREPARER}" --db-root "${db_root}")
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

if [[ "${environment}" == "hpc" && "${runtime}" == "docker" ]]; then
    die "Docker is not supported by the SLURM launcher; use --apptainer or --conda"
fi
if [[ -n "${database_config}" && ! -r "${database_config}" ]]; then
    die "database configuration is not readable: ${database_config}"
fi
if [[ "${mode}" == "run" ]]; then
    validate_forwarded_input
fi

declare -a profiles=()
case "${mode}" in
    run)
        profiles=("${environment}" "${runtime}")
        ;;
    *)
        die "unsupported execution mode '${mode}'"
        ;;
esac

profile_list="$(IFS=,; printf '%s' "${profiles[*]}")"
declare -a nextflow_command=(nextflow)
if [[ -n "${database_config}" ]]; then
    nextflow_command+=(-c "${database_config}")
fi
nextflow_command+=(run "${PIPELINE_ROOT}" -profile "${profile_list}")
if [[ "${resume}" == true ]]; then
    nextflow_command+=(-resume)
fi
nextflow_command+=("${forwarded_args[@]}")

if [[ "${dry_run}" == true ]]; then
    print_command "${nextflow_command[@]}"
    exit 0
fi

check_nextflow
check_runtime
"${nextflow_command[@]}"
