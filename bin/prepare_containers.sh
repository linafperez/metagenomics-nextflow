#!/usr/bin/env bash

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"

cache_dir=""
tmp_dir=""
validate_only=false
list_only=false
include_gpu=false

usage() {
    cat <<'USAGE'
Usage:
  ./bin/prepare_containers.sh \
      --cache-dir DIR \
      --tmp-dir DIR \
      [--validate-only] \
      [--list-only] \
      [--include-gpu]

The script:
  1. Discovers container images referenced by the pipeline.
  2. Skips valid images already present in the shared cache.
  3. Pulls missing BioContainers images from the Galaxy Singularity depot first.
  4. Falls back to docker:// when necessary.
  5. Validates every resulting SIF.
  6. Writes container_manifest.tsv with SHA-256 checksums.
USAGE
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --cache-dir)
            cache_dir="${2:?missing value for --cache-dir}"
            shift 2
            ;;
        --tmp-dir)
            tmp_dir="${2:?missing value for --tmp-dir}"
            shift 2
            ;;
        --validate-only)
            validate_only=true
            shift
            ;;
        --list-only)
            list_only=true
            shift
            ;;
        --include-gpu)
            include_gpu=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "ERROR: unknown argument: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

command -v python3 >/dev/null 2>&1 || {
    echo "ERROR: python3 was not found" >&2
    exit 2
}

mapfile -t images < <(
    python3 - "$ROOT" "$include_gpu" <<'PY'
import re
import sys
from pathlib import Path

root = Path(sys.argv[1])
include_gpu = sys.argv[2].lower() == "true"

files = []
files.extend(root.glob("modules/**/*.nf"))
files.extend(root.glob("subworkflows/**/*.nf"))
files.extend(root.glob("workflows/**/*.nf"))
files.append(root / "nextflow.config")

pattern = re.compile(
    r"((?:[A-Za-z0-9._-]+/)+[A-Za-z0-9._-]+:"
    r"[A-Za-z0-9][A-Za-z0-9._-]*|"
    r"[A-Za-z0-9._-]+:[0-9][A-Za-z0-9._-]*)"
)

images = set()

for path in files:
    if not path.is_file():
        continue

    for line in path.read_text(errors="replace").splitlines():
        if "container" not in line.lower():
            continue

        if not include_gpu and "gpu" in line.lower():
            continue

        for image in pattern.findall(line):
            image = image.removeprefix("docker://")

            # Repository-local fallback names are not the production images.
            if image.startswith("metagenomics/"):
                continue

            images.add(image)

for image in sorted(images):
    print(image)
PY
)

if [[ ${#images[@]} -eq 0 ]]; then
    echo "ERROR: no container images were discovered" >&2
    exit 2
fi

echo "Discovered ${#images[@]} container images:"
printf '  %s\n' "${images[@]}"

if [[ "$list_only" == true ]]; then
    exit 0
fi

[[ -n "$cache_dir" ]] || {
    echo "ERROR: --cache-dir is required" >&2
    exit 2
}

[[ -n "$tmp_dir" ]] || {
    echo "ERROR: --tmp-dir is required" >&2
    exit 2
}

command -v singularity >/dev/null 2>&1 || {
    echo "ERROR: singularity was not found on PATH" >&2
    exit 2
}

mkdir -p "$cache_dir" "$tmp_dir"

cache_dir="$(cd "$cache_dir" && pwd -P)"
tmp_dir="$(cd "$tmp_dir" && pwd -P)"

build_cache="${tmp_dir}/container-pull-cache.$$"
mkdir -p "$build_cache"

cleanup() {
    rm -rf "$build_cache"
}
trap cleanup EXIT

export SINGULARITY_TMPDIR="$tmp_dir"

manifest="${cache_dir}/container_manifest.tsv"
manifest_tmp="${manifest}.tmp.$$"

printf 'image\tsource\tcache_file\tsha256\tstatus\n' > "$manifest_tmp"

failures=0

cache_name() {
    local image="$1"
    image="${image#docker://}"
    image="${image//\//-}"
    image="${image//:/-}"
    printf '%s.img\n' "$image"
}

validate_image() {
    local image_file="$1"
    [[ -s "$image_file" ]] || return 1
    singularity inspect "$image_file" >/dev/null 2>&1
}

pull_source() {
    local source="$1"
    local destination="$2"
    local partial="${destination}.partial.$$"

    for attempt in 1 2 3; do
        echo "  Attempt ${attempt}: ${source}"

        rm -f "$partial"
        rm -rf "$build_cache"
        mkdir -p "$build_cache"

        export SINGULARITY_CACHEDIR="$build_cache"

        if singularity pull --force "$partial" "$source"; then
            if validate_image "$partial"; then
                mv -f "$partial" "$destination"
                return 0
            fi
        fi

        rm -f "$partial"
        sleep 20
    done

    return 1
}

for image in "${images[@]}"; do
    name="$(cache_name "$image")"
    target="${cache_dir}/${name}"
    used_source="cached"

    echo
    echo "============================================================"
    echo "Image: ${image}"
    echo "Target: ${target}"
    echo "============================================================"

    if validate_image "$target"; then
        echo "READY: valid cached image"
    else
        rm -f "$target"

        if [[ "$validate_only" == true ]]; then
            echo "MISSING: ${image}" >&2
            printf '%s\t%s\t%s\t%s\t%s\n' \
                "$image" "-" "$target" "-" "missing" >> "$manifest_tmp"
            failures=$((failures + 1))
            continue
        fi

        success=false

        if [[ "$image" == quay.io/biocontainers/* ]]; then
            biocontainer="${image#quay.io/biocontainers/}"
            galaxy_source="https://depot.galaxyproject.org/singularity/${biocontainer}"

            echo "Trying prebuilt BioContainers Singularity image first..."
            if pull_source "$galaxy_source" "$target"; then
                used_source="$galaxy_source"
                success=true
            fi
        fi

        if [[ "$success" == false ]]; then
            docker_source="docker://${image}"

            echo "Trying OCI/Docker source..."
            if pull_source "$docker_source" "$target"; then
                used_source="$docker_source"
                success=true
            fi
        fi

        if [[ "$success" == false ]]; then
            echo "ERROR: could not prepare ${image}" >&2
            printf '%s\t%s\t%s\t%s\t%s\n' \
                "$image" "-" "$target" "-" "failed" >> "$manifest_tmp"
            failures=$((failures + 1))
            continue
        fi
    fi

    sha256="$(sha256sum "$target" | awk '{print $1}')"

    printf '%s\t%s\t%s\t%s\t%s\n' \
        "$image" "$used_source" "$target" "$sha256" "ready" >> "$manifest_tmp"

    echo "READY: ${image}"
done

mv -f "$manifest_tmp" "$manifest"

echo
echo "Container manifest:"
echo "$manifest"

if (( failures > 0 )); then
    echo "ERROR: ${failures} container image(s) are not ready" >&2
    exit 1
fi

echo "ALL CONTAINERS READY"
