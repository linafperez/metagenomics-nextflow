process COMEBIN {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/comebin:1.0.4--hdfd78af_1'

    input:
    tuple val(meta), path(contigs), path(bams)

    output:
    tuple val(meta), path('*.bins/*.{fa,fna,fasta}', arity: '1..*'), emit: bins
    tuple val(meta), path('*.contigs2bin.tsv'), emit: contigs2bin
    tuple val(meta), path('*.comebin'), optional: true, emit: native_outputs
    tuple val(meta), path('*.comebin.log'), emit: log
    tuple val(meta), path('*.gpu_metrics.tsv'), optional: true, emit: gpu_metrics
    tuple val("${task.process}"), val('comebin'), val('1.0.4'), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: meta.id
    def keep_native_outputs = params.save_intermediates.toString().toBoolean()
    def bin_transfer = keep_native_outputs ? 'cp' : 'mv'
    def gpu_enabled = params.enableGpu.toString().toBoolean()
    def gpu_interval = params.gpuTelemetryInterval as int
    if (gpu_enabled && gpu_interval < 1) {
        error 'gpuTelemetryInterval must be at least one second'
    }
    def session_id = workflow.sessionId.toString()
    def attempt = task.attempt.toString()
    def gpu_metrics_file = "${prefix}.COMEBIN.__SESSION_ID__.__ATTEMPT__.gpu_metrics.tsv"
        .replace('__SESSION_ID__', session_id)
        .replace('__ATTEMPT__', attempt)
    def gpu_guard = gpu_enabled ? '''
    python3 - <<'PY'
import torch

if torch.__version__.split('+', 1)[0] != '2.1.2':
    raise SystemExit(f'COMEBin GPU mode requires PyTorch 2.1.2; found {torch.__version__}')
if torch.version.cuda != '11.8':
    raise SystemExit(f'COMEBin GPU mode requires CUDA 11.8 PyTorch; found {torch.version.cuda}')
if not torch.cuda.is_available():
    raise SystemExit('COMEBin GPU mode cannot access CUDA')
if torch.cuda.device_count() != 1:
    raise SystemExit(f'COMEBin GPU mode requires exactly one visible device; found {torch.cuda.device_count()}')
PY
    '''.stripIndent() : ''
    def gpu_monitor = gpu_enabled ? '''
    GPU_METRICS_FILE='__METRICS_FILE__'
    GPU_MONITOR_PID=''
    printf 'timestamp\tprocess\tsample_id\tsession_id\tattempt\tgpu_index\tgpu_uuid\tgpu_name\tutilization_gpu_percent\tmemory_used_mib\tmemory_total_mib\n' > "$GPU_METRICS_FILE"
    if command -v nvidia-smi >/dev/null 2>&1; then
        (
            while :; do
                timestamp="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
                nvidia-smi --query-gpu=index,uuid,name,utilization.gpu,memory.used,memory.total \
                    --format=csv,noheader,nounits 2>/dev/null \
                    | sed 's/,[[:space:]]*/,/g' \
                    | while IFS=, read -r gpu_index gpu_uuid gpu_name gpu_util gpu_memory_used gpu_memory_total; do
                        printf '%s\t__PROCESS__\t__SAMPLE_ID__\t__SESSION_ID__\t__ATTEMPT__\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                            "$timestamp" "$gpu_index" "$gpu_uuid" "$gpu_name" "$gpu_util" "$gpu_memory_used" "$gpu_memory_total"
                    done >> "$GPU_METRICS_FILE"
                sleep __INTERVAL__
            done
        ) &
        GPU_MONITOR_PID=$!
    fi
    stop_gpu_monitor() {
        if [[ -n "$GPU_MONITOR_PID" ]]; then
            kill "$GPU_MONITOR_PID" 2>/dev/null || true
            wait "$GPU_MONITOR_PID" 2>/dev/null || true
            GPU_MONITOR_PID=''
        fi
    }
    trap stop_gpu_monitor EXIT
    '''.stripIndent()
        .replace('__METRICS_FILE__', gpu_metrics_file)
        .replace('__PROCESS__', 'COMEBIN')
        .replace('__SAMPLE_ID__', meta.id.toString())
        .replace('__SESSION_ID__', session_id)
        .replace('__ATTEMPT__', attempt)
        .replace('__INTERVAL__', gpu_interval.toString()) : ''
    def gpu_finish = gpu_enabled ? 'stop_gpu_monitor\ntrap - EXIT' : ''
    def cuda_environment = gpu_enabled ? '' : "CUDA_VISIBLE_DEVICES=''"
    def bam_files = bams instanceof List ? bams : [bams]
    if (!bam_files) {
        error 'COMEBIN requires at least one reference-sorted BAM file'
    }
    def bam_links = bam_files.collect { bam -> "ln -s \"\$(readlink -f '${bam}')\" \"bam/${bam.name}\"" }.join('\n')
    def standardize_bins = '''
    shopt -s nullglob
    native_bins=("$NATIVE_BINS"/*.fa "$NATIVE_BINS"/*.fna "$NATIVE_BINS"/*.fasta)
    if [ "${#native_bins[@]}" -eq 0 ]; then
        echo "COMEBin did not produce any FASTA bins" >&2
        exit 1
    fi
    __BIN_TRANSFER__ -- "${native_bins[@]}" "$BINS_DIR/"

    : > "$MAP_FILE"
    bin_files=("$BINS_DIR"/*.fa "$BINS_DIR"/*.fna "$BINS_DIR"/*.fasta)
    for bin_file in "${bin_files[@]}"; do
        bin_name=$(basename "$bin_file")
        bin_name=${bin_name%.*}
        awk -v bin="$bin_name" 'BEGIN { OFS="\t" } /^>/ { sub(/^>/, ""); split($0, fields, /[[:space:]]+/); print fields[1], bin }' "$bin_file" >> "$MAP_FILE"
    done
    test -s "$MAP_FILE"
    '''.stripIndent().replace('__BIN_TRANSFER__', bin_transfer)

    """
    set -euo pipefail

    mkdir -p bam "${prefix}.comebin.bins"
    ${bam_links}
    ${gpu_guard}
    ${gpu_monitor}

    ${cuda_environment} run_comebin.sh \
        -a "${contigs}" \
        -o "${prefix}.comebin" \
        -p bam \
        -t ${task.cpus} \
        ${args} \
        2> >(tee "${prefix}.comebin.log" >&2)
    ${gpu_finish}

    NATIVE_BINS="${prefix}.comebin/comebin_res/comebin_res_bins"
    BINS_DIR="${prefix}.comebin.bins"
    MAP_FILE="${prefix}.comebin.contigs2bin.tsv"
    ${standardize_bins}
    rm -rf -- bam
    ${keep_native_outputs ? '' : "rm -rf -- '${prefix}.comebin'"}
    """
}
