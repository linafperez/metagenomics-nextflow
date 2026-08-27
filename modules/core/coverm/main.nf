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

}
