process GTDBTK_CLASSIFY {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${task.ext.container ?: 'quay.io/biocontainers/gtdbtk:2.6.1--pyh1f0d9b5_2'}"

    input:
    tuple val(meta), path(mags, arity: '1..*')
    path gtdbtk_db

    output:
    tuple val(meta), path('*.gtdbtk.bac120.summary.tsv'), emit: bac120_summary
    tuple val(meta), path('*.gtdbtk.ar53.summary.tsv'), emit: ar53_summary
    tuple val(meta), path('*.gtdbtk'), emit: results
    tuple val(meta), path('*.gtdbtk.log'), emit: log
    tuple val("${task.process}"), val('gtdbtk'), val('2.6.1'), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: meta.id
    def pplacer_cpus = Math.min(task.cpus as int, 8)
    def mag_files = mags instanceof List ? mags : [mags]

    if (!mag_files) {
        error 'GTDBTK_CLASSIFY requires at least one MAG FASTA file'
    }

    def entries = mag_files.collect { mag ->
        def raw_id = mag.name.replaceFirst(/(?i)\.(fa|fna|fasta)(\.gz)?$/, '')
        def mag_id = raw_id.replaceAll(/[^A-Za-z0-9_.-]/, '_')
        if (!mag_id) {
            error "GTDBTK_CLASSIFY could not derive an identifier from '${mag.name}'"
        }
        [file: mag, id: mag_id]
    }

    if (entries.collect { entry -> entry.id }.toSet().size() != entries.size()) {
        error 'GTDBTK_CLASSIFY MAG identifiers are not unique after sanitization'
    }

    def stage_commands = entries.collect { entry ->
        if (entry.file.name.toLowerCase().endsWith('.gz')) {
            return "gzip -cd \"${entry.file}\" > \"genomes/${entry.id}.fna\""
        }
        "ln -s \"\$PWD/${entry.file}\" \"genomes/${entry.id}.fna\""
    }.join('\n')

    def batch_lines = entries.collect { entry ->
        "genomes/${entry.id}.fna\t${entry.id}"
    }.join('\n')

    def summary_header = [
        'user_genome',
        'classification',
        'closest_genome_reference',
        'closest_genome_reference_radius',
        'closest_genome_taxonomy',
        'closest_genome_ani',
        'closest_genome_af',
        'closest_placement_reference',
        'closest_placement_radius',
        'closest_placement_taxonomy',
        'closest_placement_ani',
        'closest_placement_af',
        'pplacer_taxonomy',
        'classification_method',
        'note',
        'other_related_references(genome_id,species_name,radius,ANI,AF)',
        'msa_percent',
        'translation_table',
        'red_value',
        'warnings'
    ].join('\t')

    """
    set -euo pipefail

    mkdir -p genomes
    ${stage_commands}

    printf '%s\n' '${batch_lines}' > genomes.batch.tsv

    export GTDBTK_DATA_PATH="\$PWD/${gtdbtk_db}"

    echo "Running GTDB-Tk taxonomic classification"

    gtdbtk classify_wf \
        --batchfile genomes.batch.tsv \
        --out_dir "${prefix}.gtdbtk" \
        --prefix "${prefix}" \
        --cpus ${task.cpus} \
        --pplacer_cpus ${pplacer_cpus} \
        ${args} \
        2>&1 | tee "${prefix}.gtdbtk.log"

    bac_summary=\$(find "${prefix}.gtdbtk" -type f -name '*.bac120.summary.tsv' -print -quit)
    ar_summary=\$(find "${prefix}.gtdbtk" -type f -name '*.ar53.summary.tsv' -print -quit)

    if [[ -n "\${bac_summary}" ]]; then
        cp "\${bac_summary}" "${prefix}.gtdbtk.bac120.summary.tsv"
    else
        printf '%s\n' '${summary_header}' > "${prefix}.gtdbtk.bac120.summary.tsv"
    fi

    if [[ -n "\${ar_summary}" ]]; then
        cp "\${ar_summary}" "${prefix}.gtdbtk.ar53.summary.tsv"
    else
        printf '%s\n' '${summary_header}' > "${prefix}.gtdbtk.ar53.summary.tsv"
    fi
    """

    stub:
    def prefix    = task.ext.prefix ?: meta.id
    def mag_files = mags instanceof List ? mags : [mags]
    def entries   = mag_files.collect { mag ->
        def raw_id = mag.name.replaceFirst(/(?i)\.(fa|fna|fasta)(\.gz)?$/, '')
        [id: raw_id.replaceAll(/[^A-Za-z0-9_.-]/, '_')]
    }
    def summary_header = [
        'user_genome', 'classification', 'closest_genome_reference', 'closest_genome_reference_radius',
        'closest_genome_taxonomy', 'closest_genome_ani', 'closest_genome_af', 'closest_placement_reference',
        'closest_placement_radius', 'closest_placement_taxonomy', 'closest_placement_ani',
        'closest_placement_af', 'pplacer_taxonomy', 'classification_method', 'note',
        'other_related_references(genome_id,species_name,radius,ANI,AF)', 'msa_percent',
        'translation_table', 'red_value', 'warnings'
    ].join('\t')
    def bacterial_rows = entries.collect { entry ->
        [
            entry.id,
            'd__Bacteria;p__Stubphylum;c__Stubclass;o__Stuborder;f__Stubfamily;g__Stubgenus;s__',
            'N/A', 'N/A', 'N/A', 'N/A', 'N/A',
            'N/A', 'N/A', 'N/A', 'N/A', 'N/A',
            'd__Bacteria;p__Stubphylum;c__Stubclass;o__Stuborder;f__Stubfamily;g__Stubgenus;s__',
            'taxonomic classification fully defined by topology',
            'N/A', 'N/A', '100.0', '11', '0.5', 'Stub classification'
        ].join('\t')
    }.join('\n')

    """
    set -euo pipefail

    mkdir -p "${prefix}.gtdbtk/classify" "${prefix}.gtdbtk/identify" "${prefix}.gtdbtk/align"

    printf '%s\n' '${summary_header}' > "${prefix}.gtdbtk/${prefix}.bac120.summary.tsv"
    printf '%s\n' '${bacterial_rows}' >> "${prefix}.gtdbtk/${prefix}.bac120.summary.tsv"

    printf '%s\n' '${summary_header}' > "${prefix}.gtdbtk/${prefix}.ar53.summary.tsv"

    cp "${prefix}.gtdbtk/${prefix}.bac120.summary.tsv" "${prefix}.gtdbtk.bac120.summary.tsv"
    cp "${prefix}.gtdbtk/${prefix}.ar53.summary.tsv" "${prefix}.gtdbtk.ar53.summary.tsv"

    printf '%s\n' '{"version":"2.6.1","database_version":"r226","status":"stub"}' \
        > "${prefix}.gtdbtk/gtdbtk.json"

    printf 'GTDB-Tk stub classification completed\n' > "${prefix}.gtdbtk/gtdbtk.log"
    printf 'GTDB-Tk stub classification completed\n' > "${prefix}.gtdbtk.log"
    """
}
