process COVERM_CONTIG {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/coverm:0.7.0--hcb7b614_4'

    input:
    tuple val(meta), path(contigs), val(sample_ids), path(reads)

    output:
    tuple val(meta), path('*.metabat_depth.tsv'), emit: metabat_depth
    tuple val(meta), path('*.vamb_abundance.tsv'), emit: vamb_abundance
    tuple val(meta), path('bam/*.bam'), emit: bams
    tuple val(meta), path('*.coverm.log'), emit: log
    tuple val("${task.process}"), val('coverm'), val('0.7.0'), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def coverage_args = task.ext.coverage_args ?: ''
    def prefix        = task.ext.prefix ?: meta.id
    def ids           = sample_ids instanceof List ? sample_ids : [sample_ids]
    def read_files    = reads instanceof List ? reads.flatten() : [reads]

    if (read_files.size() != ids.size() * 2) {
        error "COVERM_CONTIG expected two reads for each of ${ids.size()} sample IDs, but received ${read_files.size()} files"
    }

    def safe_ids = ids.collect { id -> id.toString().replaceAll(/[^A-Za-z0-9_.-]/, '_') }
    if (safe_ids.toSet().size() != safe_ids.size()) {
        error 'COVERM_CONTIG sample IDs are not unique after filename sanitization'
    }

    def pairs = safe_ids.withIndex().collect { id, index ->
        def r1 = read_files[index * 2]
        def r2 = read_files[index * 2 + 1]
        def r1_ext = r1.name.toLowerCase().endsWith('.gz') ? '.fastq.gz' : '.fastq'
        def r2_ext = r2.name.toLowerCase().endsWith('.gz') ? '.fastq.gz' : '.fastq'
        [id: id, r1: r1, r2: r2, r1_link: "reads/${id}_R1${r1_ext}", r2_link: "reads/${id}_R2${r2_ext}"]
    }

    def link_commands = pairs.collect { pair ->
        "ln -s \"\$(readlink -f '${pair.r1}')\" \"${pair.r1_link}\"\nln -s \"\$(readlink -f '${pair.r2}')\" \"${pair.r2_link}\""
    }.join('\n')
    def coupled_args = pairs.collect { pair -> "\"${pair.r1_link}\" \"${pair.r2_link}\"" }.join(' ')
    def rename_bams = pairs.collect { pair ->
        [
            "cache_matches=(cache/*\"${pair.id}_R1\"*.bam)",
            'if [ "${#cache_matches[@]}" -ne 1 ]; then',
            "    echo \"Expected one cached BAM for sample ${pair.id}, found \${#cache_matches[@]}\" >&2",
            '    exit 1',
            'fi',
            "mv -- \"\${cache_matches[0]}\" \"bam/${pair.id}.bam\""
        ].join('\n')
    }.join('\n')
    def bam_args = safe_ids.collect { id -> "\"bam/${id}.bam\"" }.join(' ')
    def vamb_header = (['contigname'] + safe_ids).join('\\t')
    def expected_columns = safe_ids.size() + 1

    """
    set -euo pipefail

    mkdir -p reads cache bam
    ${link_commands}

    coverm contig \
        --coupled ${coupled_args} \
        --reference "${contigs}" \
        --mapper minimap2-sr \
        --methods mean \
        --threads ${task.cpus} \
        --output-format dense \
        --output-file "${prefix}.mapping_mean.tsv" \
        --bam-file-cache-directory cache \
        ${args} \
        2> >(tee "${prefix}.coverm.log" >&2)

    shopt -s nullglob
    ${rename_bams}

    coverm contig \
        --bam-files ${bam_args} \
        --methods mean \
        --threads ${task.cpus} \
        --output-format dense \
        --output-file "${prefix}.mean_depth.raw.tsv" \
        ${coverage_args} \
        2> >(tee -a "${prefix}.coverm.log" >&2)

    awk -v expected=${expected_columns} -v header="${vamb_header}" '
        BEGIN { FS=OFS="\t" }
        NR == 1 {
            if (NF != expected) {
                printf "Unexpected CoverM mean-depth column count: %d (expected %d)\n", NF, expected > "/dev/stderr"
                exit 1
            }
            print header
            next
        }
        { print }
    ' "${prefix}.mean_depth.raw.tsv" > "${prefix}.vamb_abundance.tsv"

    coverm contig \
        --bam-files ${bam_args} \
        --methods metabat \
        --threads ${task.cpus} \
        --output-format dense \
        --output-file "${prefix}.metabat_depth.tsv" \
        ${coverage_args} \
        2> >(tee -a "${prefix}.coverm.log" >&2)

    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    def ids = sample_ids instanceof List ? sample_ids : [sample_ids]
    def safe_ids = ids.collect { id -> id.toString().replaceAll(/[^A-Za-z0-9_.-]/, '_') }
    def bam_args = safe_ids.collect { id -> "\"${id}\"" }.join(' ')
    def vamb_header = (['contigname'] + safe_ids).join('\\t')
    def vamb_values = safe_ids.collect { '0' }.join('\\t')
    def metabat_header = (['contigName', 'contigLen', 'totalAvgDepth'] + safe_ids.collectMany { id -> ["${id}.bam", "${id}.bam-var"] }).join('\\t')
    def metabat_values = (['1000', '0'] + safe_ids.collectMany { ['0', '0'] }).join('\\t')

    """
    set -euo pipefail

    mkdir -p bam
    python3 - ${bam_args} <<'PY'
    import struct
    import sys
    import zlib
    from pathlib import Path

    sam_header = b"@HD\\tVN:1.6\\tSO:coordinate\\n@SQ\\tSN:stub_contig\\tLN:1000\\n"
    reference_name = b"stub_contig\\0"
    raw_bam = (
        b"BAM\\1"
        + struct.pack("<i", len(sam_header))
        + sam_header
        + struct.pack("<i", 1)
        + struct.pack("<i", len(reference_name))
        + reference_name
        + struct.pack("<i", 1000)
    )
    compressor = zlib.compressobj(level=6, wbits=-15)
    payload = compressor.compress(raw_bam) + compressor.flush()
    block_size = 18 + len(payload) + 8
    bgzf_header = bytes.fromhex("1f8b08040000000000ff060042430200") + struct.pack(
        "<H", block_size - 1
    )
    trailer = struct.pack("<II", zlib.crc32(raw_bam) & 0xFFFFFFFF, len(raw_bam))
    eof_block = bytes.fromhex(
        "1f8b08040000000000ff0600424302001b0003000000000000000000"
    )

    for sample_id in sys.argv[1:]:
        Path("bam", f"{sample_id}.bam").write_bytes(
            bgzf_header + payload + trailer + eof_block
        )
    PY
    contig_id=\$(awk '/^>/{sub(/^>/, ""); split(\$0, fields, /[[:space:]]+/); print fields[1]; exit}' "${contigs}")
    contig_id=\${contig_id:-stub_contig}

    printf '%b\n' "${vamb_header}" > "${prefix}.vamb_abundance.tsv"
    printf '%s\t%b\n' "\${contig_id}" "${vamb_values}" >> "${prefix}.vamb_abundance.tsv"
    printf '%b\n' "${metabat_header}" > "${prefix}.metabat_depth.tsv"
    printf '%s\t%b\n' "\${contig_id}" "${metabat_values}" >> "${prefix}.metabat_depth.tsv"
    printf 'CoverM contig stub completed\n' > "${prefix}.coverm.log"

    """
}


process COVERM_GENOME {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/coverm:0.7.0--hcb7b614_4'

    input:
    tuple val(meta), path(mags), val(sample_ids), path(reads)

    output:
    tuple val(meta), path('*.mag_abundance.tsv'), emit: abundance
    tuple val(meta), path('*.coverm.log'), emit: log
    tuple val("${task.process}"), val('coverm'), val('0.7.0'), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: meta.id
    def ids        = sample_ids instanceof List ? sample_ids : [sample_ids]
    def read_files = reads instanceof List ? reads.flatten() : [reads]
    def mag_files  = mags instanceof List ? mags : [mags]

    if (read_files.size() != ids.size() * 2) {
        error "COVERM_GENOME expected two reads for each of ${ids.size()} sample IDs, but received ${read_files.size()} files"
    }
    if (!mag_files) {
        error 'COVERM_GENOME requires at least one MAG FASTA file'
    }

    def safe_ids = ids.collect { id -> id.toString().replaceAll(/[^A-Za-z0-9_.-]/, '_') }
    if (safe_ids.toSet().size() != safe_ids.size()) {
        error 'COVERM_GENOME sample IDs are not unique after filename sanitization'
    }

    def pairs = safe_ids.withIndex().collect { id, index ->
        def r1 = read_files[index * 2]
        def r2 = read_files[index * 2 + 1]
        def r1_ext = r1.name.toLowerCase().endsWith('.gz') ? '.fastq.gz' : '.fastq'
        def r2_ext = r2.name.toLowerCase().endsWith('.gz') ? '.fastq.gz' : '.fastq'
        [id: id, r1: r1, r2: r2, r1_link: "reads/${id}_R1${r1_ext}", r2_link: "reads/${id}_R2${r2_ext}"]
    }
    def link_commands = pairs.collect { pair ->
        "ln -s \"\$(readlink -f '${pair.r1}')\" \"${pair.r1_link}\"\nln -s \"\$(readlink -f '${pair.r2}')\" \"${pair.r2_link}\""
    }.join('\n')
    def coupled_args = pairs.collect { pair -> "\"${pair.r1_link}\" \"${pair.r2_link}\"" }.join(' ')
    def mag_args = mag_files.collect { mag -> "\"${mag}\"" }.join(' ')

    """
    set -euo pipefail

    mkdir -p reads
    ${link_commands}

    coverm genome \
        --genome-fasta-files ${mag_args} \
        --coupled ${coupled_args} \
        --mapper minimap2-sr \
        --methods relative_abundance mean covered_bases length \
        --proper-pairs-only \
        --min-read-percent-identity 95 \
        --min-read-aligned-percent 75 \
        --min-covered-fraction 0 \
        --contig-end-exclusion 0 \
        --threads ${task.cpus} \
        --output-format dense \
        --output-file "${prefix}.mag_abundance.raw.tsv" \
        ${args} \
        2> >(tee "${prefix}.coverm.log" >&2)

    awk '
        BEGIN { FS=OFS="\t" }
        NR == 1 {
            for (column = 2; column <= NF; column++) {
                sub(/^.*\\//, "", \$column)
                sub(/_R1\\.fastq(\\.gz)? /, " ", \$column)
                sub(/ Covered Bases\$/, " Covered Fraction", \$column)
            }
            print
            next
        }
        {
            for (column = 4; column <= NF; column += 4) {
                length_column = column + 1
                \$column = \$length_column > 0 ? \$column / \$length_column : 0
            }
        }
        { print }
    ' "${prefix}.mag_abundance.raw.tsv" > "${prefix}.mag_abundance.tsv"
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    def ids = sample_ids instanceof List ? sample_ids : [sample_ids]
    def safe_ids = ids.collect { id -> id.toString().replaceAll(/[^A-Za-z0-9_.-]/, '_') }
    def mag_files = mags instanceof List ? mags : [mags]
    def header = (['Genome'] + safe_ids.collectMany { id -> ["${id} Relative Abundance (%)", "${id} Mean", "${id} Covered Fraction", "${id} Length"] }).join('\\t')
    def values = safe_ids.collectMany { ['0', '0', '0', '0'] }.join('\\t')
    def rows = mag_files.collect { mag ->
        def mag_id = mag.name.replaceFirst(/\.[^.]+$/, '')
        "printf '%s\\t%b\\n' \"${mag_id}\" \"${values}\" >> \"${prefix}.mag_abundance.tsv\""
    }.join('\n')

    """
    set -euo pipefail

    printf '%b\n' "${header}" > "${prefix}.mag_abundance.tsv"
    ${rows}
    printf 'CoverM genome stub completed\n' > "${prefix}.coverm.log"

    """
}
