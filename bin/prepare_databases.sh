#!/usr/bin/env bash

set -euo pipefail
umask 0027

readonly PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
readonly GENCODE_RELEASE="44"
readonly GENCODE_ASSEMBLY="GRCh38.p14"
readonly GENCODE_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.p14.genome.fa.gz"
readonly GTDBTK_RELEASE="226"
readonly GTDBTK_URL="https://data.gtdb.ecogenomic.org/releases/release226/226.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r226_data.tar.gz"
readonly GTDBTK_MD5_URL="https://data.gtdb.ecogenomic.org/releases/release226/MD5SUM.txt"
readonly PHYLOPHLAN_URL="https://cmprod1.cibio.unitn.it/databases/PhyloPhlAn/phylophlan.tar"
readonly PHYLOPHLAN_MD5_URL="https://cmprod1.cibio.unitn.it/databases/PhyloPhlAn/phylophlan.md5"
readonly INTERPROSCAN_VERSION="5.59-91.0"
readonly INTERPROSCAN_URL="https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.59-91.0/alt/interproscan-data-5.59-91.0.tar.gz"
readonly INTERPROSCAN_MD5_URL="https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.59-91.0/alt/interproscan-data-5.59-91.0.tar.gz.md5"
readonly CHECKM2_VERSION="1.1.0"
readonly GUNC_VERSION="1.0.6"
readonly EGGNOG_MAPPER_VERSION="2.1.13"
readonly EGGNOG_DATABASE_VERSION="5.0.2"
readonly EGGNOG_BASE_URL="https://downloads.eggnogdb.org/emapper/emapperdb-5.0.2"
readonly PHYLOPHLAN_VERSION="3.1.1"

readonly -a RESOURCE_NAMES=(
    human_reference
    bowtie2_index
    checkm2
    gunc
    gtdbtk
    phylophlan
    eggnog
    interproscan
    genemarks2
)

db_root=""
human_reference_source=""
genemark_home=""
genemark_key=""
bowtie2_build_command="bowtie2-build"
checkm2_command="checkm2"
gunc_command="gunc"
eggnog_downloader=""
jobs=8
dry_run=false
check_only=false
redownload_all=false
selection_restricted=false
current_stage=""

declare -A selected=()
declare -A skipped=()
declare -A redownload=()
declare -A status=()
declare -A detail=()
declare -A resource_path=()
declare -A provenance=()
declare -A version=()

usage() {
    cat <<'USAGE'
Usage:
  prepare_databases.sh --db-root PATH [options]

Required:
  --db-root PATH             Shared database root outside the Git repository.

Selection and validation:
  --only NAME[,NAME...]      Prepare only the named resources; repeatable.
  --skip NAME[,NAME...]      Leave named resources unchanged; repeatable.
  --redownload NAME[,NAME...] Rebuild named resources with recoverable backups.
  --redownload-all           Rebuild every selected resource.
  --check-only               Validate resources and regenerate metadata only.
  --dry-run                  Print planned actions without downloads or writes.
  --jobs N                   Threads for Bowtie2 index construction (default: 8).

External or site-provided inputs:
  --human-reference FILE     Use a local GRCh38.p14 FASTA or FASTA.GZ.
  --genemark-home DIR        Licensed GeneMarkS-2 1.15 installation.
  --genemark-key FILE        Licensed GeneMark key; it is never copied.

Tool command overrides:
  --bowtie2-build PATH       Bowtie2 index builder (default: bowtie2-build).
  --checkm2-command PATH     CheckM2 executable (default: checkm2).
  --gunc-command PATH        GUNC executable (default: gunc).
  --eggnog-downloader PATH   Optional patched eggNOG 5.0.2 downloader override.

Resource names:
  human_reference, bowtie2_index, checkm2, gunc, gtdbtk, phylophlan,
  eggnog, interproscan, genemarks2

The script resumes HTTP transfers, validates expected files, preserves replaced
resources as timestamped backups, and writes database_manifest.tsv plus
metagenomics_databases.config under --db-root. GTDB-Tk and InterProScan are very
large. GeneMarkS-2 is never downloaded because its license requires the user to
obtain both the software and key directly from GeneMark.
USAGE
}

log() {
    printf '%s\n' "$*"
}

warn() {
    printf 'Warning: %s\n' "$*" >&2
}

die() {
    printf 'Error: %s\n' "$*" >&2
    exit 2
}

require_value() {
    local option="$1"
    local value="${2:-}"
    [[ -n "${value}" && "${value}" != --* ]] || die "${option} requires a value"
}

normalize_resource_name() {
    local name="${1,,}"
    name="${name//-/_}"
    case "${name}" in
        human|reference|human_reference) printf 'human_reference\n' ;;
        bowtie2|host_index|bowtie2_index) printf 'bowtie2_index\n' ;;
        checkm2|gunc|gtdbtk|phylophlan|eggnog|interproscan|genemarks2)
            printf '%s\n' "${name}"
            ;;
        genemark|genemark_s2) printf 'genemarks2\n' ;;
        *) die "unknown resource name '${1}'" ;;
    esac
}

add_resource_list() {
    local target_name="$1"
    local raw_list="$2"
    local raw_name normalized
    local -a names=()

    IFS=, read -r -a names <<< "${raw_list}"
    for raw_name in "${names[@]}"; do
        [[ -n "${raw_name}" ]] || die "empty resource name in '${raw_list}'"
        normalized="$(normalize_resource_name "${raw_name}")"
        case "${target_name}" in
            selected) selected["${normalized}"]=true ;;
            skipped) skipped["${normalized}"]=true ;;
            redownload) redownload["${normalized}"]=true ;;
            *) die "internal resource-list error" ;;
        esac
    done
}

canonical_path() {
    realpath -m -- "$1"
}

path_is_within() {
    local child="$1"
    local parent="$2"
    [[ "${child}" == "${parent}" || "${child}" == "${parent}/"* ]]
}

validate_db_root() {
    local project_canonical db_canonical probe git_root

    command -v realpath >/dev/null 2>&1 || die "realpath is required"
    db_canonical="$(canonical_path "${db_root}")"
    project_canonical="$(canonical_path "${PROJECT_ROOT}")"
    [[ "${db_canonical}" != "/" ]] || die "--db-root cannot be the filesystem root"
    [[ "${db_canonical}" != "/tmp" ]] || die "--db-root must be a dedicated directory"
    path_is_within "${db_canonical}" "${project_canonical}" \
        && die "--db-root must be outside the pipeline Git repository"

    probe="${db_canonical}"
    while [[ ! -e "${probe}" && "${probe}" != "/" ]]; do
        probe="$(dirname "${probe}")"
    done
    git_root="$({ git -C "${probe}" rev-parse --show-toplevel 2>/dev/null || true; })"
    if [[ -n "${git_root}" ]]; then
        git_root="$(canonical_path "${git_root}")"
        path_is_within "${db_canonical}" "${git_root}" \
            && die "--db-root must not be inside a Git worktree"
    fi
    db_root="${db_canonical}"
}

command_available() {
    local executable="$1"
    if [[ "${executable}" == */* ]]; then
        [[ -x "${executable}" ]]
    else
        command -v "${executable}" >/dev/null 2>&1
    fi
}

print_command() {
    local argument
    printf 'Command:'
    for argument in "$@"; do
        printf ' %q' "${argument}"
    done
    printf '\n'
}

run_command() {
    if [[ "${dry_run}" == true ]]; then
        print_command "$@"
        return 0
    fi
    "$@"
}

cleanup_stage() {
    if [[ -n "${current_stage}" && -e "${current_stage}" ]]; then
        case "${current_stage}" in
            "${db_root}"/.staging/*) rm -rf -- "${current_stage}" ;;
            *) warn "refusing to clean unexpected staging path ${current_stage}" ;;
        esac
    fi
}

trap cleanup_stage EXIT INT TERM

new_stage() {
    local resource="$1"
    local stage_parent="${db_root}/.staging"

    if [[ "${dry_run}" == true ]]; then
        current_stage="${stage_parent}/${resource}.DRY_RUN"
        printf '%s\n' "${current_stage}"
        return
    fi
    mkdir -p -- "${stage_parent}"
    current_stage="$(mktemp -d "${stage_parent}/${resource}.XXXXXXXX")"
    chmod 0750 "${current_stage}"
    printf '%s\n' "${current_stage}"
}

backup_existing() {
    local target="$1"
    local timestamp backup

    [[ -e "${target}" || -L "${target}" ]] || return 0
    timestamp="$(date -u +'%Y%m%dT%H%M%SZ')"
    backup="${target}.backup.${timestamp}"
    if [[ "${dry_run}" == true ]]; then
        print_command mv -- "${target}" "${backup}"
    else
        mv -- "${target}" "${backup}"
        log "Preserved previous resource at ${backup}"
    fi
}

install_staged_directory() {
    local staged="$1"
    local target="$2"

    case "${staged}" in
        "${db_root}"/.staging/*) ;;
        *) die "refusing to install an unexpected staging path: ${staged}" ;;
    esac
    backup_existing "${target}"
    if [[ "${dry_run}" == true ]]; then
        print_command mv -- "${staged}" "${target}"
    else
        mkdir -p -- "$(dirname "${target}")"
        mv -- "${staged}" "${target}"
    fi
    current_stage=""
}

download_file() {
    local url="$1"
    local output="$2"

    if [[ "${dry_run}" == true ]]; then
        print_command curl --fail --location --continue-at - --retry 5 \
            --retry-all-errors --output "${output}" "${url}"
        return
    fi
    command -v curl >/dev/null 2>&1 || die "curl is required for database downloads"
    mkdir -p -- "$(dirname "${output}")"
    curl --fail --location --continue-at - --retry 5 --retry-all-errors \
        --output "${output}" "${url}"
}

prepare_cache_for_redownload() {
    local resource="$1"
    local cache_file="$2"
    local timestamp

    if [[ "${redownload_all}" != true && "${redownload[${resource}]:-false}" != true ]]; then
        return
    fi
    [[ -e "${cache_file}" ]] || return
    timestamp="$(date -u +'%Y%m%dT%H%M%SZ')"
    run_command mv -- "${cache_file}" "${cache_file}.previous.${timestamp}"
}

verify_md5_file() {
    local archive="$1"
    local checksum_file="$2"
    local archive_name expected actual

    if [[ "${dry_run}" == true ]]; then
        log "Would verify MD5 for ${archive}"
        return
    fi
    command -v md5sum >/dev/null 2>&1 || die "md5sum is required"
    archive_name="$(basename "${archive}")"
    expected="$(grep -F "${archive_name}" "${checksum_file}" | head -n 1 | awk '{print $1}')"
    if [[ -z "${expected}" ]]; then
        expected="$(awk 'NF {print $1; exit}' "${checksum_file}")"
    fi
    [[ "${expected}" =~ ^[0-9a-fA-F]{32}$ ]] \
        || die "no valid MD5 checksum was found for ${archive_name}"
    actual="$(md5sum "${archive}" | awk '{print $1}')"
    [[ "${actual,,}" == "${expected,,}" ]] \
        || die "MD5 validation failed for ${archive_name}"
}

validate_fasta() {
    local fasta="$1"
    [[ -s "${fasta}" ]] || return 1
    awk 'BEGIN {valid=0} /^>/ {valid=1; exit} END {exit valid ? 0 : 1}' "${fasta}"
}

validate_bowtie2_index() {
    local prefix="$1"
    local extension suffix
    for extension in bt2 bt2l; do
        for suffix in 1 2 3 4 rev.1 rev.2; do
            [[ -s "${prefix}.${suffix}.${extension}" ]] || break
        done
        [[ "${suffix}" == "rev.2" && -s "${prefix}.rev.2.${extension}" ]] && return 0
    done
    return 1
}

find_checkm2_database() {
    find -L "${resource_path[checkm2]}" -type f -name 'uniref100.KO.1.dmnd' -print -quit 2>/dev/null
}

find_gunc_database() {
    find -L "${resource_path[gunc]}" -type f -name '*.dmnd' -print -quit 2>/dev/null
}

validate_checkm2() {
    local database_file
    database_file="$(find_checkm2_database)"
    [[ -n "${database_file}" && -s "${database_file}" ]]
}

validate_gunc() {
    local database_file
    database_file="$(find_gunc_database)"
    [[ -n "${database_file}" && -s "${database_file}" ]]
}

validate_gtdbtk() {
    local directory="${resource_path[gtdbtk]}"
    [[ -d "${directory}/markers" \
        && -d "${directory}/metadata" \
        && -d "${directory}/masks" ]]
}

validate_phylophlan() {
    local database_file
    [[ -d "${resource_path[phylophlan]}" ]] || return 1
    database_file="$(find -L "${resource_path[phylophlan]}" -maxdepth 2 \
        -type f \( -name 'phylophlan*.faa' -o -name 'phylophlan*.dmnd' \) \
        -print -quit 2>/dev/null)"
    [[ -n "${database_file}" && -s "${database_file}" ]]
}

validate_eggnog() {
    [[ -s "${resource_path[eggnog]}/eggnog.db" \
        && -s "${resource_path[eggnog]}/eggnog_proteins.dmnd" \
        && -s "${resource_path[eggnog]}/eggnog.taxa.db" ]]
}

validate_interproscan() {
    local entry
    [[ -d "${resource_path[interproscan]}" ]] || return 1
    entry="$(find -L "${resource_path[interproscan]}" -mindepth 1 -maxdepth 2 \
        -type f -size +0c -print -quit 2>/dev/null)"
    [[ -n "${entry}" ]]
}

validate_genemarks2() {
    [[ -n "${genemark_home}" && -d "${genemark_home}" \
        && -n "${genemark_key}" && -s "${genemark_key}" ]] || return 1
    find -L "${genemark_home}" -type f -name 'gms2.pl' -print -quit 2>/dev/null \
        | grep -q .
}

resource_valid() {
    case "$1" in
        human_reference) validate_fasta "${resource_path[human_reference]}" ;;
        bowtie2_index) validate_bowtie2_index "${resource_path[bowtie2_index]}" ;;
        checkm2) validate_checkm2 ;;
        gunc) validate_gunc ;;
        gtdbtk) validate_gtdbtk ;;
        phylophlan) validate_phylophlan ;;
        eggnog) validate_eggnog ;;
        interproscan) validate_interproscan ;;
        genemarks2) validate_genemarks2 ;;
        *) return 1 ;;
    esac
}

should_prepare() {
    local resource="$1"

    [[ "${selected[${resource}]:-false}" == true ]] || return 1
    [[ "${skipped[${resource}]:-false}" != true ]] || return 1
    [[ "${check_only}" != true ]] || return 1
    if [[ "${redownload_all}" == true || "${redownload[${resource}]:-false}" == true ]]; then
        return 0
    fi
    ! resource_valid "${resource}"
}

prepare_human_reference() {
    local stage archive reference_output
    if ! should_prepare human_reference; then
        return 0
    fi
    log "Preparing GRCh38.p14 human reference"
    new_stage human_reference >/dev/null
    stage="${current_stage}"
    reference_output="${stage}/GRCh38.p14.fa"

    if [[ -n "${human_reference_source}" ]]; then
        [[ -s "${human_reference_source}" ]] \
            || die "human reference was not found: ${human_reference_source}"
        if [[ "${human_reference_source}" == *.gz ]]; then
            run_command bash -o pipefail -c 'gzip -cd -- "$1" > "$2"' \
                prepare-reference "${human_reference_source}" "${reference_output}"
        else
            run_command cp -- "${human_reference_source}" "${reference_output}"
        fi
    else
        archive="${db_root}/.downloads/GRCh38.p14.genome.fa.gz"
        prepare_cache_for_redownload human_reference "${archive}"
        download_file "${GENCODE_URL}" "${archive}"
        run_command gzip -t "${archive}"
        run_command bash -o pipefail -c 'gzip -cd -- "$1" > "$2"' \
            prepare-reference "${archive}" "${reference_output}"
    fi
    if [[ "${dry_run}" != true ]]; then
        validate_fasta "${reference_output}" || die "GRCh38.p14 FASTA validation failed"
    fi
    install_staged_directory "${stage}" "$(dirname "${resource_path[human_reference]}")"
}

prepare_bowtie2_index() {
    local stage prefix
    if ! should_prepare bowtie2_index; then
        return 0
    fi
    if [[ "${dry_run}" != true ]]; then
        resource_valid human_reference \
            || die "Bowtie2 index preparation requires a valid GRCh38.p14 reference"
        command_available "${bowtie2_build_command}" \
            || die "Bowtie2 index builder is unavailable: ${bowtie2_build_command}"
    fi
    log "Building GRCh38.p14 Bowtie2 index"
    new_stage bowtie2_index >/dev/null
    stage="${current_stage}"
    prefix="${stage}/GRCh38_p14"
    run_command "${bowtie2_build_command}" --threads "${jobs}" \
        "${resource_path[human_reference]}" "${prefix}"
    if [[ "${dry_run}" != true ]]; then
        validate_bowtie2_index "${prefix}" || die "Bowtie2 index validation failed"
    fi
    install_staged_directory "${stage}" "$(dirname "${resource_path[bowtie2_index]}")"
}

prepare_checkm2() {
    local stage staged_database
    if ! should_prepare checkm2; then
        return 0
    fi
    if [[ "${dry_run}" != true ]]; then
        command_available "${checkm2_command}" \
            || die "CheckM2 is unavailable; install CheckM2 ${CHECKM2_VERSION} or use --checkm2-command"
    fi
    log "Downloading the CheckM2 ${CHECKM2_VERSION} database"
    new_stage checkm2 >/dev/null
    stage="${current_stage}"
    run_command "${checkm2_command}" database --download --path "${stage}"
    if [[ "${dry_run}" != true ]]; then
        staged_database="$(find -L "${stage}" -type f -name 'uniref100.KO.1.dmnd' -print -quit)"
        [[ -n "${staged_database}" && -s "${staged_database}" ]] \
            || die "CheckM2 database validation failed"
    fi
    install_staged_directory "${stage}" "${resource_path[checkm2]}"
}

prepare_gunc() {
    local stage staged_database
    if ! should_prepare gunc; then
        return 0
    fi
    if [[ "${dry_run}" != true ]]; then
        command_available "${gunc_command}" \
            || die "GUNC is unavailable; install GUNC ${GUNC_VERSION} or use --gunc-command"
    fi
    log "Downloading the GUNC ProGenomes 2.1 database"
    new_stage gunc >/dev/null
    stage="${current_stage}"
    run_command "${gunc_command}" download_db --database progenomes "${stage}"
    if [[ "${dry_run}" != true ]]; then
        staged_database="$(find -L "${stage}" -type f -name '*.dmnd' -print -quit)"
        [[ -n "${staged_database}" && -s "${staged_database}" ]] \
            || die "GUNC database validation failed"
    fi
    install_staged_directory "${stage}" "${resource_path[gunc]}"
}

prepare_gtdbtk() {
    local stage archive checksum extracted candidate
    if ! should_prepare gtdbtk; then
        return 0
    fi
    log "Downloading GTDB-Tk release ${GTDBTK_RELEASE}"
    archive="${db_root}/.downloads/gtdbtk_r226_data.tar.gz"
    checksum="${db_root}/.downloads/GTDB_r226_MD5SUM.txt"
    prepare_cache_for_redownload gtdbtk "${archive}"
    download_file "${GTDBTK_URL}" "${archive}"
    download_file "${GTDBTK_MD5_URL}" "${checksum}"
    verify_md5_file "${archive}" "${checksum}"
    new_stage gtdbtk >/dev/null
    stage="${current_stage}"
    run_command tar -xzf "${archive}" -C "${stage}"

    extracted="${stage}"
    if [[ "${dry_run}" != true && ! -d "${extracted}/markers" ]]; then
        candidate="$(find "${stage}" -mindepth 1 -maxdepth 3 -type d -name markers \
            -printf '%h\n' -quit)"
        [[ -n "${candidate}" ]] || die "GTDB-Tk release ${GTDBTK_RELEASE} validation failed"
        extracted="${candidate}"
    fi
    install_staged_directory "${extracted}" "${resource_path[gtdbtk]}"
    if [[ "${dry_run}" != true && "${extracted}" != "${stage}" ]]; then
        current_stage="${stage}"
        cleanup_stage
        current_stage=""
    fi
}

prepare_phylophlan() {
    local stage archive checksum database_file extracted
    if ! should_prepare phylophlan; then
        return 0
    fi
    log "Downloading the official PhyloPhlAn marker database"
    archive="${db_root}/.downloads/phylophlan.tar"
    checksum="${db_root}/.downloads/phylophlan.md5"
    prepare_cache_for_redownload phylophlan "${archive}"
    download_file "${PHYLOPHLAN_URL}" "${archive}"
    download_file "${PHYLOPHLAN_MD5_URL}" "${checksum}"
    verify_md5_file "${archive}" "${checksum}"
    new_stage phylophlan >/dev/null
    stage="${current_stage}"
    run_command tar -xf "${archive}" -C "${stage}"
    extracted="${stage}/phylophlan"
    if [[ "${dry_run}" != true ]]; then
        database_file="$(find "${stage}" -maxdepth 3 -type f \
            -name 'phylophlan*.faa' -size +0c -print -quit)"
        [[ -n "${database_file}" ]] || die "PhyloPhlAn marker database validation failed"
        extracted="$(dirname "${database_file}")"
    fi
    install_staged_directory "${extracted}" "${resource_path[phylophlan]}"
    if [[ "${dry_run}" != true && "${extracted}" != "${stage}" ]]; then
        current_stage="${stage}"
        cleanup_stage
        current_stage=""
    fi
}

prepare_eggnog() {
    local stage annotation_archive taxonomy_archive diamond_archive
    if ! should_prepare eggnog; then
        return 0
    fi
    if [[ -n "${eggnog_downloader}" && "${dry_run}" != true ]]; then
        command_available "${eggnog_downloader}" \
            || die "configured eggNOG downloader is unavailable: ${eggnog_downloader}"
    fi
    log "Downloading the eggNOG ${EGGNOG_DATABASE_VERSION} mapper database"
    new_stage eggnog >/dev/null
    stage="${current_stage}"
    if [[ -n "${eggnog_downloader}" ]]; then
        run_command "${eggnog_downloader}" -y --data_dir "${stage}"
    else
        annotation_archive="${db_root}/.downloads/eggnog-${EGGNOG_DATABASE_VERSION}.db.gz"
        taxonomy_archive="${db_root}/.downloads/eggnog-${EGGNOG_DATABASE_VERSION}.taxa.tar.gz"
        diamond_archive="${db_root}/.downloads/eggnog-${EGGNOG_DATABASE_VERSION}.proteins.dmnd.gz"
        prepare_cache_for_redownload eggnog "${annotation_archive}"
        prepare_cache_for_redownload eggnog "${taxonomy_archive}"
        prepare_cache_for_redownload eggnog "${diamond_archive}"
        download_file "${EGGNOG_BASE_URL}/eggnog.db.gz" "${annotation_archive}"
        download_file "${EGGNOG_BASE_URL}/eggnog.taxa.tar.gz" "${taxonomy_archive}"
        download_file "${EGGNOG_BASE_URL}/eggnog_proteins.dmnd.gz" "${diamond_archive}"
        run_command cp -- "${annotation_archive}" "${stage}/eggnog.db.gz"
        run_command cp -- "${diamond_archive}" "${stage}/eggnog_proteins.dmnd.gz"
        run_command gzip --decompress --force "${stage}/eggnog.db.gz"
        run_command gzip --decompress --force "${stage}/eggnog_proteins.dmnd.gz"
        run_command tar -xzf "${taxonomy_archive}" -C "${stage}"
    fi
    if [[ "${dry_run}" != true ]]; then
        [[ -s "${stage}/eggnog.db" \
            && -s "${stage}/eggnog_proteins.dmnd" \
            && -s "${stage}/eggnog.taxa.db" ]] \
            || die "eggNOG database validation failed"
    fi
    install_staged_directory "${stage}" "${resource_path[eggnog]}"
}

prepare_interproscan() {
    local stage archive checksum data_directory final_directory
    if ! should_prepare interproscan; then
        return 0
    fi
    log "Downloading InterProScan ${INTERPROSCAN_VERSION} data"
    archive="${db_root}/.downloads/interproscan-data-${INTERPROSCAN_VERSION}.tar.gz"
    checksum="${db_root}/.downloads/interproscan-data-${INTERPROSCAN_VERSION}.tar.gz.md5"
    prepare_cache_for_redownload interproscan "${archive}"
    download_file "${INTERPROSCAN_URL}" "${archive}"
    download_file "${INTERPROSCAN_MD5_URL}" "${checksum}"
    verify_md5_file "${archive}" "${checksum}"
    new_stage interproscan >/dev/null
    stage="${current_stage}"
    run_command tar -xzf "${archive}" -C "${stage}"
    final_directory="${stage}/prepared_data"
    if [[ "${dry_run}" == true ]]; then
        print_command mv -- "${stage}/interproscan-${INTERPROSCAN_VERSION}/data" "${final_directory}"
    else
        data_directory="$(find "${stage}" -mindepth 2 -maxdepth 4 -type d -name data -print -quit)"
        [[ -n "${data_directory}" ]] || die "InterProScan data directory was not found"
        find "${data_directory}" -mindepth 1 -maxdepth 2 -type f -size +0c -print -quit \
            | grep -q . || die "InterProScan data validation failed"
        mv -- "${data_directory}" "${final_directory}"
    fi
    install_staged_directory "${final_directory}" "${resource_path[interproscan]}"
    if [[ "${dry_run}" != true ]]; then
        current_stage="${stage}"
        cleanup_stage
        current_stage=""
    fi
}

inspect_resources() {
    local resource
    for resource in "${RESOURCE_NAMES[@]}"; do
        if [[ "${selected[${resource}]:-false}" != true ]]; then
            if resource_valid "${resource}"; then
                status["${resource}"]="ready"
                detail["${resource}"]="Valid existing resource; not selected for preparation"
            else
                status["${resource}"]="not_selected"
                detail["${resource}"]="Not requested"
            fi
        elif [[ "${skipped[${resource}]:-false}" == true ]]; then
            if resource_valid "${resource}"; then
                status["${resource}"]="ready"
                detail["${resource}"]="Valid existing resource; preparation skipped"
            else
                status["${resource}"]="skipped"
                detail["${resource}"]="Preparation explicitly skipped"
            fi
        elif resource_valid "${resource}"; then
            status["${resource}"]="ready"
            detail["${resource}"]="Validated"
        elif [[ "${resource}" == "genemarks2" ]]; then
            status["${resource}"]="license_required"
            detail["${resource}"]="Provide --genemark-home and --genemark-key"
        else
            status["${resource}"]="missing"
            detail["${resource}"]="Expected files were not found"
        fi
    done
}

groovy_value() {
    local value="$1"
    if [[ -z "${value}" ]]; then
        printf 'null'
        return
    fi
    value="${value//\\/\\\\}"
    value="${value//\'/\\\'}"
    printf "'%s'" "${value}"
}

config_path_if_ready() {
    local resource="$1"
    local value="$2"
    if [[ "${status[${resource}]}" == "ready" ]]; then
        printf '%s' "${value}"
    fi
}

write_metadata() {
    local checked_at manifest config temporary_manifest temporary_config
    local resource gunc_file checkm2_path host_index gtdb_path phylo_path eggnog_path interpro_path
    local configured_genemark_home configured_genemark_key

    checked_at="$(date -u +'%Y-%m-%dT%H:%M:%SZ')"
    manifest="${db_root}/database_manifest.tsv"
    config="${db_root}/metagenomics_databases.config"
    if [[ "${dry_run}" == true ]]; then
        log "Would write ${manifest}"
        log "Would write ${config}"
        return
    fi

    temporary_manifest="${db_root}/.database_manifest.tsv.$$"
    temporary_config="${db_root}/.metagenomics_databases.config.$$"
    printf 'resource\tversion\tstatus\tpath\tprovenance\tchecked_at\tdetail\n' \
        > "${temporary_manifest}"
    for resource in "${RESOURCE_NAMES[@]}"; do
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "${resource}" "${version[${resource}]}" "${status[${resource}]}" \
            "${resource_path[${resource}]}" "${provenance[${resource}]}" \
            "${checked_at}" "${detail[${resource}]}" >> "${temporary_manifest}"
    done

    host_index="$(config_path_if_ready bowtie2_index "${resource_path[bowtie2_index]}")"
    checkm2_path="$(config_path_if_ready checkm2 "${resource_path[checkm2]}")"
    gunc_file=""
    if [[ "${status[gunc]}" == "ready" ]]; then
        gunc_file="$(find_gunc_database)"
    fi
    gtdb_path="$(config_path_if_ready gtdbtk "${resource_path[gtdbtk]}")"
    phylo_path="$(config_path_if_ready phylophlan "${resource_path[phylophlan]}")"
    eggnog_path="$(config_path_if_ready eggnog "${resource_path[eggnog]}")"
    interpro_path="$(config_path_if_ready interproscan "${resource_path[interproscan]}")"
    configured_genemark_home="$(config_path_if_ready genemarks2 "${genemark_home}")"
    configured_genemark_key="$(config_path_if_ready genemarks2 "${genemark_key}")"

    {
        printf '// Generated by bin/prepare_databases.sh at %s\n' "${checked_at}"
        printf 'params {\n'
        printf '    host_bowtie2_index = %s\n' "$(groovy_value "${host_index}")"
        printf '    checkm2_db         = %s\n' "$(groovy_value "${checkm2_path}")"
        printf '    gunc_db            = %s\n' "$(groovy_value "${gunc_file}")"
        printf '    gtdbtk_db          = %s\n' "$(groovy_value "${gtdb_path}")"
        printf '    phylophlan_db      = %s\n' "$(groovy_value "${phylo_path}")"
        printf '    eggnog_db          = %s\n' "$(groovy_value "${eggnog_path}")"
        printf '    interproscan_data  = %s\n' "$(groovy_value "${interpro_path}")"
        printf '    genemark_home      = %s\n' "$(groovy_value "${configured_genemark_home}")"
        printf '    genemark_key       = %s\n' "$(groovy_value "${configured_genemark_key}")"
        printf '}\n'
    } > "${temporary_config}"

    chmod 0640 "${temporary_manifest}" "${temporary_config}"
    mv -f -- "${temporary_manifest}" "${manifest}"
    mv -f -- "${temporary_config}" "${config}"
    log "Wrote ${manifest}"
    log "Wrote ${config}"
}

while (($#)); do
    case "$1" in
        --db-root)
            require_value "$1" "${2:-}"
            db_root="$2"
            shift
            ;;
        --db-root=*) db_root="${1#*=}" ;;
        --only)
            require_value "$1" "${2:-}"
            selection_restricted=true
            add_resource_list selected "$2"
            shift
            ;;
        --only=*)
            selection_restricted=true
            add_resource_list selected "${1#*=}"
            ;;
        --skip)
            require_value "$1" "${2:-}"
            add_resource_list skipped "$2"
            shift
            ;;
        --skip=*) add_resource_list skipped "${1#*=}" ;;
        --redownload)
            require_value "$1" "${2:-}"
            add_resource_list redownload "$2"
            shift
            ;;
        --redownload=*) add_resource_list redownload "${1#*=}" ;;
        --redownload-all) redownload_all=true ;;
        --check-only) check_only=true ;;
        --dry-run) dry_run=true ;;
        --jobs)
            require_value "$1" "${2:-}"
            jobs="$2"
            shift
            ;;
        --jobs=*) jobs="${1#*=}" ;;
        --human-reference)
            require_value "$1" "${2:-}"
            human_reference_source="$2"
            shift
            ;;
        --human-reference=*) human_reference_source="${1#*=}" ;;
        --genemark-home)
            require_value "$1" "${2:-}"
            genemark_home="$2"
            shift
            ;;
        --genemark-home=*) genemark_home="${1#*=}" ;;
        --genemark-key)
            require_value "$1" "${2:-}"
            genemark_key="$2"
            shift
            ;;
        --genemark-key=*) genemark_key="${1#*=}" ;;
        --bowtie2-build)
            require_value "$1" "${2:-}"
            bowtie2_build_command="$2"
            shift
            ;;
        --bowtie2-build=*) bowtie2_build_command="${1#*=}" ;;
        --checkm2-command)
            require_value "$1" "${2:-}"
            checkm2_command="$2"
            shift
            ;;
        --checkm2-command=*) checkm2_command="${1#*=}" ;;
        --gunc-command)
            require_value "$1" "${2:-}"
            gunc_command="$2"
            shift
            ;;
        --gunc-command=*) gunc_command="${1#*=}" ;;
        --eggnog-downloader)
            require_value "$1" "${2:-}"
            eggnog_downloader="$2"
            shift
            ;;
        --eggnog-downloader=*) eggnog_downloader="${1#*=}" ;;
        -h|--help)
            usage
            exit 0
            ;;
        *) die "unknown option '$1'" ;;
    esac
    shift
done

[[ -n "${db_root}" ]] || die "--db-root PATH is required"
[[ "${jobs}" =~ ^[1-9][0-9]*$ ]] || die "--jobs must be a positive integer"
if [[ "${check_only}" == true \
    && ("${redownload_all}" == true || ${#redownload[@]} -gt 0) ]]; then
    die "--check-only cannot be combined with redownload options"
fi
validate_db_root

if [[ "${selection_restricted}" == false ]]; then
    for resource in "${RESOURCE_NAMES[@]}"; do
        selected["${resource}"]=true
    done
fi

for resource in "${!redownload[@]}"; do
    [[ "${selected[${resource}]:-false}" == true ]] \
        || die "--redownload resource '${resource}' is not selected"
    [[ "${skipped[${resource}]:-false}" != true ]] \
        || die "resource '${resource}' cannot be both skipped and redownloaded"
done

if [[ -n "${human_reference_source}" ]]; then
    human_reference_source="$(canonical_path "${human_reference_source}")"
fi
if [[ -n "${genemark_home}" ]]; then
    genemark_home="$(canonical_path "${genemark_home}")"
fi
if [[ -n "${genemark_key}" ]]; then
    genemark_key="$(canonical_path "${genemark_key}")"
fi

resource_path[human_reference]="${db_root}/human/GRCh38.p14/GRCh38.p14.fa"
resource_path[bowtie2_index]="${db_root}/human/bowtie2/GRCh38_p14"
resource_path[checkm2]="${db_root}/checkm2/1.1.0"
resource_path[gunc]="${db_root}/gunc/1.0.6/progenomes_2.1"
resource_path[gtdbtk]="${db_root}/gtdbtk/release226"
resource_path[phylophlan]="${db_root}/phylophlan/phylophlan"
resource_path[eggnog]="${db_root}/eggnog/5.0.2"
resource_path[interproscan]="${db_root}/interproscan/5.59-91.0/data"
resource_path[genemarks2]="${genemark_home:-external_license_required}"

version[human_reference]="GENCODE release ${GENCODE_RELEASE}; ${GENCODE_ASSEMBLY}"
version[bowtie2_index]="Bowtie2 2.5.4 index; ${GENCODE_ASSEMBLY}"
version[checkm2]="${CHECKM2_VERSION}"
version[gunc]="${GUNC_VERSION}; ProGenomes 2.1"
version[gtdbtk]="release ${GTDBTK_RELEASE}"
version[phylophlan]="${PHYLOPHLAN_VERSION} marker database"
version[eggnog]="${EGGNOG_DATABASE_VERSION}; mapper ${EGGNOG_MAPPER_VERSION}"
version[interproscan]="${INTERPROSCAN_VERSION}"
version[genemarks2]="1.15"

provenance[human_reference]="${GENCODE_URL}"
provenance[bowtie2_index]="Locally built from the GENCODE release ${GENCODE_RELEASE} reference"
provenance[checkm2]="checkm2 database --download --path"
provenance[gunc]="gunc download_db --database progenomes; ProGenomes 2.1 database"
provenance[gtdbtk]="${GTDBTK_URL}"
provenance[phylophlan]="${PHYLOPHLAN_URL}"
provenance[eggnog]="${EGGNOG_BASE_URL}; eggNOG-mapper ${EGGNOG_MAPPER_VERSION} compatible database"
provenance[interproscan]="${INTERPROSCAN_URL}"
provenance[genemarks2]="User-provided licensed installation and key"

if [[ -n "${human_reference_source}" ]]; then
    provenance[human_reference]="User-provided GRCh38.p14 FASTA: ${human_reference_source}"
fi

if [[ "${dry_run}" != true ]]; then
    mkdir -p -- "${db_root}/.downloads" "${db_root}/.staging"
    if command -v flock >/dev/null 2>&1; then
        exec 9>"${db_root}/.prepare_databases.lock"
        flock -n 9 || die "another database preparation process is using ${db_root}"
    fi
fi

prepare_human_reference
prepare_bowtie2_index
prepare_checkm2
prepare_gunc
prepare_gtdbtk
prepare_phylophlan
prepare_eggnog
prepare_interproscan

if [[ "${selected[genemarks2]:-false}" == true \
    && "${skipped[genemarks2]:-false}" != true \
    && ! validate_genemarks2 ]]; then
    warn "GeneMarkS-2 remains pending; obtain version 1.15 and its key from GeneMark"
fi

inspect_resources
write_metadata

missing_count=0
for resource in "${RESOURCE_NAMES[@]}"; do
    printf '%-18s %s\n' "${resource}" "${status[${resource}]}"
    if [[ "${selected[${resource}]:-false}" == true \
        && "${skipped[${resource}]:-false}" != true \
        && "${status[${resource}]}" != "ready" ]]; then
        ((missing_count += 1))
    fi
done

if [[ "${check_only}" == true && "${dry_run}" != true && ${missing_count} -gt 0 ]]; then
    die "${missing_count} selected resource(s) failed validation"
fi

if [[ "${missing_count}" -gt 0 ]]; then
    warn "${missing_count} selected resource(s) remain pending; see database_manifest.tsv"
fi
